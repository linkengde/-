# 六个液态 Ag pilot：GPAW runner 离线复用审查

2026-10-11，Asia/Shanghai。只阅读公开源码、规则及安装包源码；没有构造 GPAW calculator，没有运行 MPI、DFT、MD，也没有修改 A/B 文件。当前用户已授权六个构型生成、六个 DFT pilot 和必要的同结构 k 网格对照；其余 24 项仍 HOLD。运行前须满足冻结执行包中的来源、实际坐标、QA、持久工作盘及资源门槛。

## 可复用来源

仓库根 `/tmp/cloud-c-v19-design`，以下均为相对路径：

- `research/mace-v12-transfer/coordination/reports/window-c/unified_potential_dft_data_plan_20261010/METHOD_AND_ROLE_RULES.md`：两种能量口径、PAW/PBE/PW500、SCF、数据角色规则。
- `research/mace-v12-transfer/coordination/reports/window-c/five_cloud_dft_parallel_plan_20261010/HANDOFF_D.md`：液体 Gamma/2x2x2 数值预检和两项 D 候选。
- `research/mace-v12-transfer/coordination/reports/window-c/five_cloud_preflight_and_liquid_pilot_readiness_20261011/DISK_AND_RESOURCE_BUDGET.md`：真实历史大 checkpoint 证据及 28.8–29.0 GB 风险。
- `research/mace-v12-transfer/coordination/reports/window-b/v19_pure_phase_runner_and_bulk_mesh_preparation/run_pw_reference.py`：新 GPAW API hook、MPI 先初始化后 ASE 读、4 ranks、输出完整力与17有效数字、原子顺序、收敛钩子。
- 同目录 `entry.py`：输入方法/hash/order/PBC 校验、MPI rank binding 下 CPU quota 判断。
- 同目录 `verify_result.py`：完整结果清单、几何/PBC/order/ID精确保持、有限 Nx3 力、log/summary/progress 三重 SCF 核验。
- 同目录 `run_queue.py`：独占锁及 one-job 模式的思路；不能复用自动 `sync_tasks.claim/progress/publish` 调用，因为那会登记或修改已有 A/B 任务。

本机 GPAW26.7源码只读路径：

- `/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw/convergence_criteria.py` 第122/179/218行；
- `.../gpaw/new/ase_interface.py` 第171–195、285–294、336–360行；
- `.../gpaw/new/gpw.py` 第93–136行；
- `.../gpaw/new/builder.py` 第518–567行；
- `.../gpaw/new/calculation.py` 第277–285行。

## 科学设置和 SCF 的准确语义

冻结完整配置为 GPAW26.7.0、ASE3.29.0、PAW data1.2.1、Ag.PBE.gz已归档hash、PBE、`PW(500)`、`FermiDirac(0.10)`、`convergence={'energy':1e-5,'density':1e-5,'eigenstates':1e-5}`、`maxiter=160`。旧 runner 的 `Mixer(0.05,8,100)` 和 `parallel={'sl_auto':True}` 可以作为一致性候选，但也须写入执行包方法身份而不是隐藏默认值。液体 PBC=[True,True,True]；无 slab dipole solver，`poissonsolver={}`。不增加几何优化、外场、MD步骤、stress标签要求、charge/spin变化。

三项 1e-5 不是同一个单位：

- energy：默认最近三步 **extrapolated total energy 的 peak-to-peak / valence electron**，单位 eV/valence-electron。
- density：积分绝对密度误差 / valence electron，单位 electrons/valence-electron。
- eigenstates：Kohn–Sham 残差平方积分 / valence electron，单位 eV²/valence-electron。

因此不能写成“总能量收敛 1e-5 eV/体系”。继续保留原参数，不在新任务中替换成 energy-per-atom或不同绝对阈值。实际日志须显示三项及其 GPAW 描述、收敛步数和 converged marker。

`from gpaw import GPAW` 在26.7调用新 API；旧源码 `gpaw.calculator.GPAW` 是 old API。新 runner 应禁止转到 old API，以确保 `calc.hooks['scf_step']`/`['converged']`语义一致。

## 能量与力输出

- `atoms.get_potential_energy(force_consistent=False)` 返回 `energy`，即 native extrapolated energy。
- `atoms.get_potential_energy(force_consistent=True)` 返回 `free_energy`，即与有限展宽原子力一致的能量。
- 保留 `energy_eV_cell`、`free_energy_eV_cell`、`energy_key='GPAW energy / force_consistent=False'`、`force_consistent_energy_key='GPAW free_energy / force_consistent=True'` 和 `forces_eV_A` 完整 Nx3 数组，不能只保留 fmax/RMS。
- free_energy必须存在且有限；不要沿用旧runner允许 `None` 的弱条件。
- 输入和结果保持相同 atom ID/order、positions、cell、PBC，不包裹坐标、不重排、不静默改变几何。输出17有效数字，可用 JSON 向量+extxyz；extxyz显式给 free_energy 和 native_energy，不让ASE SinglePoint自动属性覆盖命名。
- 新六个标签仍是候选数据，不能自动变成训练数据；现有15 Ag/10 TSC纯相控制角色不变。公开已评分父族不是独立测试。

## k 点定义与数值 gate

冻结 `kpts=(1,1,1)` 为真正 Gamma；`kpts=(2,2,2)` 是 standard Monkhorst–Pack shifted even mesh，不要在不同云端混用 gamma-centered2x2x2或不记录offset。归档实际 BZ/IBZ 坐标、权重、time-reversal/symmetry配置及 NIBZ。random液体结构通常没有点群；time reversal 可能减少IBZ但不能只凭推测用于预算。Gamma实dtype与2meshcomplexdtype也记录。

在**同一精确几何 SHA、同一原子顺序**上比较 Gamma vs2mesh；不要比较不同seed/密度几何的能量或力。所有非k点参数完全相同：

- `dE_native_meV_per_atom = 1000*abs(E_native_high-E_native_low)/N`；
- `dE_free_meV_per_atom = 1000*abs(FE_high-FE_low)/N`；
- `force_vector_RMSE = sqrt(mean_i(sum_xyz((F_high-F_low)^2)))`，不是每分量均值；
- `max_atom_force_delta = max_i(norm(F_high[i]-F_low[i]))`。

数值控制通过须两种能量差各<=2 meV/atom、向量RMSE<=0.01 eV/Å、最大原子力差<=0.03 eV/Å，均SCF收敛且完整验收。定义符号/比较方向也在JSON写清。该gate检验数值离散误差，不能证明采样势精度或凝固性能。仅一个snapshot通过不能推广到三种密度全部；至少每个密度一种代表snapshot先检查；如果策略只批准一个，应保守将其他密度仍标数值待审。额外同结构mesh runs是数值控制，不新增父族或构型数量。

若2mesh不通过，停止本密度后续发布；3mesh/更密必须先停止、重新估算容量／耗时并取得用户专项批准；当前执行包只准Gamma与2mesh，不能用Gamma省资源强行给标签。若未来专项批准的更密对照通过，再由C审定可接受标签；保留旧结果为数值控制，当前不运行3mesh。

## 带数、内存与磁盘

Ag.PBE.gz源码XML的Ag是17价电子。64 Ag、neutral、无自旋时1088电子。GPAW26.7默认bands公式 `ceil(1.2*(nvalence+M)/2)+4` 后取setup cap，因此在M=0且setup bound-state cap不限制时 **657 bands**，而不是旧示意图简单翻倍的662。Ag的bound-state cap为12/atom→768，故657是本源码/已检查setup下的预测默认；仍须从本次实际初始化日志再次核验nbands、ncomponents、dtype、NIBZ、Npw、FFT和projector维度。不要自行降低band数。发生默认配置、setup、charge、spin、band padding不同则停止重新预算。

首阶波函数存储为 `itemsize*Nbands*Ncomponents*NIBZ*Npw`；GPAW内存额外包含projectors、densities、Hamiltonian、SCALAPACKwork和writer/gather缓冲；checkpoint bytes不能当RAM峰值。旧48atomTSC checkpoint7.36 GB说明“结果CSV很小”不能说明工作盘足够。

使用冻结密度/盒/mesh的实际初始化维度出预算，不使用未核实9.3 g/cm³生成盒。正式启动前必须记录cgroup `memory.max-memory.current`与MemAvailable取有效较小者；每host最多一个4rank，OMP/BLAS/MKL线程各1，无oversubscribe。内存>=16GiB是screening，不是3mesh无限放行依据。保留内存安全余量并监控memory.peak/oom.events；OOM后不自动重启。

磁盘每个task独立 numerical_run_id 和 checkpoint ledger，记录剩余bytes/inodes/quota，retained finals+current+previous+temporary三份峰值+日志输入+reserve。原30GB建议保留；低于30GB须由C依据实际输入／checkpoint／retained files／持久性／quota与reserve审定有限pilot预算，不能由runner默认放宽。30GB是建议，不改建议也不据此禁止小型采样；缺容量／持久证据时HOLD，不假称PASS，也不重复请求已授权的整轮计算许可。容量预算需要覆盖六个标签及每host对应mesh控制累积，不可从A/B清checkpoint腾空间。

## 原子checkpoint政策（旧runner必须修正）

旧runner直接 `calc.write(state.gpw,mode='all')` 每20SCF：无atomic替换、无previous、仅legacy固定reserve，没有12h限时；这些是不可继承的风险。

新政策建议：

1. 每次SCF safe hook由rank0核查实际空间/预算/elapsed，状态通过MPI broadcast让全部rank采取一致分支。不可只有rank0 raise而其他rank等待collective。
2. 所有rank同时 `calc.write(unique_tmp_path,mode='all',precision='double')`，不能只rank0调用，因为写wavefunction涉及collective通信，源码末尾有barrier。
3. GPAW写返回后rank0 fsync文件，检查ULMheader/方法/几何身份（只读metadata，不重跑），计算文件SHA和字节数；filesystem同卷保证rename原子。失败保持原current/previous及failed temp，不继续。
4. 先前已有current只有在新temp完整返回、fsync、SHA/checkpoint manifest齐备后才可成为previous，然后新temp replace current，再fsync目录。三份峰值预留；过程中哪怕中断始终有已完成generation及匹配SHA记录。所有rank barrier后继续。
5. 两代rotation若需要移除更旧的checkpoint，必须在冻结执行包中明确允许仅**新pilot自己的**已核验checkpoint retention policy；用户原话禁止清理时不要自动移除。A/B旧checkpoint绝不进入此政策。
6. convergence后forces/energies验证完再全rank写final checkpoint，compact结果另行原子归档。失败/未收敛archive不可标PASS；上次合法checkpoint只支持显式审核恢复，不支持自动restart。

SHA证明文件身份，不单独证明SCF收敛；恢复前必须核对task+inputSHA+geometry/methodSHA+PAW hashes+mesh+rank/env和ULM完整性。各mesh不共享checkpoint。不得读取A/B运行目录或清理其他host输出。

## 超时、SCF与queue行为

- `maxiter=160`是SCF硬上限，达到而未converged需捕获失败并归档FAIL/PRESERVED；不能当成有标签。
- 第一项运行elapsed>8h：允许本次作业在12h上限内继续到收敛/安全停点，但**暂停后续队列启动**等待预算审查；写 HOLD marker，queue每下一项必须读它。
- 12h在SCF safe hook检查并停止，保留最后完整checkpoint；不能声称hook是硬wallclock（单步或MPI会卡）。外部watchdog保证进程时间上限；尽量先标stop_requested在safe hook停，若rank卡死再TERM process-group，仍保存last_good而非直接覆写。不要trap中强制calc.write，因为MPI可能只部分rank响应导致死锁。KILL只作宽限结束，temp只标不完整，禁止resume该temp。
- 独占锁 `${DFT_WORK_ROOT}/.pilot.lock`，任务与run_id目录`exist_ok=False`。已有目录/进程时停止，不自动cleanup/restart。
- C可收集提交结果，但runner不自动claim或修改原`sync_tasks`A/B卡、不自动推main或force push。各host只在自己的branch提交小型公开结果+manifest+SHA，不把大.gpw提交Git。
- 明确失败原因：未知结构/无合法来源、重复任务/input、错元素数/顺序/PBC、input/method/setup SHA变化、资源/预算不足、NaN/Inf、SCF不收敛、超时、writer异常、checkpoint/持久路径不一致。失效一项并暂停后续，不静默改cutoff/kmesh/mixer/smearing/bands/geometry。
