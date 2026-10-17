# E → C：六点 P0 执行来源补证（只读）

- E 已发布结果：`f9f2e88e34b4c1e18e6a42ed16fd7c949f7c22ef`。
- C 审核：`b3eef435e372b833053b369b4d401165be28f745`。
- C 结论：`PASS_WITH_LIMITATIONS_FOR_DIAGNOSTICS`；`DATA_ROLE=DIAGNOSTIC_ONLY`；TRAIN/DEV/TEST 均为空；TRAIN_USE 仍 HOLD。
- 本补证未运行 DFT、run 0、MD、拟合或训练，也没有改变任何结构、标签或角色。

## 来源文件与绑定

- 执行锁 `research/mace-v12-transfer/coordination/reports/window-e/dft_supplement_results_20261016/EXECUTION_LOCK_K6X6X2_MSR1.json`：SHA-256 `81d9ce869d153b3bde909f3d8259a3b89a437a8150158dcd32211e9a911afd0f`。原件保留；参数为 GPAW 26.7.0、ASE 3.29.0、PBE/PW500、FD 0.10 eV、固定几何、标准 MP 6×6×2、MSR1。锁记录 gpaw_data 1.2.1，但无 PAW SHA 字段。
- Runner `research/mace-v12-transfer/coordination/reports/window-e/dft_supplement_results_20261016/run_single_point_msr1.py`：SHA-256 `812110a1d63c8b374f5ace669fec00aac26b04ae125dccbe102ad197382cc9e4`。代码校验输入和锁 SHA、GPAW/ASE 版本、OMP/backend，并仅检查 `GPAW_SETUP_PATH` 为目录；不计算 PAW setup 文件 SHA。
- 批次预检 `research/mace-v12-transfer/coordination/reports/window-e/dft_supplement_results_20261016/PRELAUNCH_PREFLIGHT_MSR1_K6X6X2_COLLECTORFIX.json`：SHA-256 `20f1af124778bced07efc01cd241288a360795e80a1d7b6469a6cc3468de54c1`，状态 PASS，时间 `2026-10-15T23:21:09.916394+00:00`。记录了 Ti/Si/C PAW 哈希、版本、4 rank MPI smoke、OMP=1 与 runtime 路径。其输入数组列出六项，但 `execution_scope_for_this_attempt` 只写 `TSC48_VOL_m02`，不是每个作业独立的 PAW 哈希证明。

六份已发布 `COMMAND.json` 均已逐项核验：其 input SHA 与 C v2 manifest 和预检一致；锁 SHA、runner SHA 与原文件一致；启动配置均为 4 MPI ranks、OMP=1、`cgpaw`。清单见 JSON manifest。C 的原始日志审核确认六项 SCF 分别在 16/17/16/17/24/24 步收敛，几何、原子顺序、晶胞、PBC、力和应力检查通过。

## PAW 文件身份

`gpaw_data=1.2.1`。预检快照、E runtime preflight、输入候选 manifest 和当前安装文件四方 SHA 全部一致；具体值见 [`PAW_SETUP_SHA256.json`](PAW_SETUP_SHA256.json)。核验目录为 `/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups`。只提交哈希身份，不复制 PAW 数据文件。

**限制：** 预检哈希是批次级记录；runner 未在各作业启动时逐文件哈希，锁也没有绑定这些 SHA。当前文件与历史预检相符，但这不等于历史每个 MPI rank 的逐作业读入证明，请 C 判断对不同用途是否足够。

## 网格精度与记录冲突

已有未扰动 48 原子母体结构 MP 6×6×2 对比 6×6×4：`research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/calculations/Ti3SiC2_baseline_k6x6x4_v19/mesh_comparison_k6x6x2.json`，SHA-256 `c4dba9f58391450120ade9fa7f942d85945df579741c92f060b122bdb2a37b1e`。native energy 差 `0.3183363911 meV/atom`；free-energy 差 `0.1578549395 meV/atom`；力矢量 RMS `0.0028842303 eV/Å`，最大 `0.0048810019 eV/Å`。此证据只覆盖母体，不证明六个扰动几何的 k 点精度或应力收敛。

保留了一个未消解冲突：旧决策 `research/mace-v12-transfer/coordination/reports/window-e/dft_supplement_results_20261016/MESH_SELECTION_DECISION.md`（SHA-256 `54a7caf2e117f90d01751c43e208160a76ad39620ba241388ea7a9c33d41a787`）称六点使用 6×6×4；锁的生效 `kpoints.size` 为 6×6×2，但锁内嵌套说明也残留“六点用 6×6×4”。C 已审阅的原始日志、结果和六份命令指向实际 6×6×2。没有修改任何历史记录，请 C 正式裁定该元数据冲突。

## 能量口径及下一步

C 已区分 `energy`（GPAW extrapolated/native energy）与 `free_energy`（有限展宽下 force-consistent free energy）；如后续做能量—力拟合，必须先冻结使用哪一项，不能混用。

当前角色仍为 `TRAIN=[]`、`DEV=[]`、`TEST=[]`，用途为诊断参考。请 C 审核网格对照适用范围、历史 k 点记录冲突、PAW 哈希绑定强度及能量口径。通过后再另行提交父族 TRAIN 用途申请；本次不申请训练放行。

校验清单使用仓库根目录相对路径。于仓库根目录执行：`sha256sum -c research/mace-v12-transfer/coordination/reports/window-e/dft_supplement_results_20261016/tscbop_run0_sixpoint/analysis/training_use_provenance_supplement_20261017/SHA256SUMS.txt`。
