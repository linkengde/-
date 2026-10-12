# R2 SCF 根因诊断与单点恢复申请

状态：`FAILED_SCF`。本报告只读分析 R1/R2 已归档结构与 GPAW 日志，并整理已做的三次 R2 诊断；没有启动计算，也没有修改原 Phase2 文件。

## 结论

当前最有证据的障碍是 R2 的电子态/密度自洽反馈不稳定，尤其是本征态残差长期停留在高位且总能反复大幅跳变。三种 Pulay mixer 设置都没有形成稳定下降轨迹。几何检查未发现 R2 原子重叠或明显异常短键，R2 的最近界面 Ag–Si 距离反而比成功收敛的 R1 更长；已记录的停止状态也不支持墙钟、内存或 MPI 资源故障是 SCF 失败原因。

不能仅凭这些日志断言电子不稳定的微观原因。R2 未收敛，所以没有最终费米能级、占据数或可恢复的电子态 checkpoint；目前最直接的区分实验是增加每次 SCF 中的平面波本征求解工作量，同时保持已锁定的 mixer 与全部物理设置不变。

## R1 与 R2 的结构及输入对比

| 项目 | R1（成功） | R2（失败） | 结论 |
|---|---|---|---|
| 原子数与组成 | 20；Ag7 C4 Si3 Ti6 | 20；Ag7 C4 Si3 Ti6 | 相同 |
| 晶胞（Å） | `[[3.071,0,0],[-1.5355,2.6595640150,0],[0,0,71.2682660002]]` | 相同 | 相同 |
| PBC | `[true,true,false]` | 相同 | 相同 |
| 真空 | 底部 18.000 Å；顶部 18.000 Å | 相同 | 相同 |
| k 点 | 8×8×1，Monkhorst–Pack 默认约定 | 相同 | 相同 |
| XC / PAW / cutoff | PBE / 同一 PAW / 500 eV | 相同 | 相同 |
| smearing / 收敛阈值 | Fermi–Dirac 0.10 eV；energy、density、eigenstates 均 1e-5 | 相同 | 相同 |
| Mixer / eigensolver | Mixer β=0.05、nmaxold=8、weight=100；PPCG niter=2 | 相同 | 相同 |
| 晶体来源 | COD-9009647 / Kisi 1998 来源族 | 相同 | 同一来源族 |
| registry | top / OT-like | bridge / MT-like | 唯一结构差异是 Ag 面内 registry 位移 |
| 最近全体系原子对 | Ti–C，2.094781 Å | Ti–C，2.094781 Å | 相同；无 R2 特有短原子对 |
| 最近 Ag–基底原子对 | Ag–Si，2.600000 Å | Ag–Si，3.019563 Å | R2 不是更短的界面接触 |
| 垂直界面间距 | 2.600000 Å | 2.600000 Å | 相同 |

两份结构中 Ti₃SiC₂ 基底坐标相同；Ag 层的面内位置按 registry 平移，原子 z 坐标不变。`R1/input.py` 与 `R2/input.py` 的文本差异只有 registry ID、来源文件路径和绑定的结构 SHA。两份运行日志均为 1 个对称操作、64 个 BZ 点/32 个 IBZ 点、136 个 bands、219 个价电子和 `[16,16,375]` 实空间网格。因此，未发现 k 点、晶胞、PBC、组成或其他已锁定电子参数的意外差异。

输入文件 SHA：R1 `d3f7b5542f071cbc432ade3277a22b0af5277cfd539d1eb48a0d006e7250ca6b`；R2 `7c1327f403c94240e6181da990fe7b078f420cc32f92a68bee3282792b0bcc66`。结构 SHA：R1 `b758b39571a6ce8b24f687e9b341c4328336c108d761b4adaa25e3fb51932890`；R2 `cb2e16376887f7cdeb05c7ddf4aacfacacbea0a000837a3888a0f3ea91785f70`。参数锁 SHA 为 `d97de3fe5a26ab6529d5e6e8ceef15adc1f3a2f548f03b59f9f092361c2c0c78`。

## SCF 残差与停止原因

`R2_SCF_RESIDUAL_COMPARISON.csv` 保存 R1、Phase2 原始 R2 和三份恢复副本的全部逐迭代能量、本征态残差、密度残差及原始日志 SHA。CSV 中的迭代能量仅用于观察收敛轨迹，**不是 DFT 标签**。

- **R1**：152 步收敛。初期能量也曾大幅波动，但后期稳定到 −188.921492 eV；最终日志将 energy、本征态和 density 判据都标记为收敛。本征态残差末值为 −7.60，density 残差末值为 −5.00（保留 GPAW 表格原始值）。费米能级 `2.004 eV` 和占据表只在收敛摘要中出现。
- **Phase2 原始 R2**：运行到配置上限 160 步后抛出 `gpaw.KohnShamConvergenceError`。最后一行为能量 `−1742.016847 eV`、eigenstate 残差 `+1.80`、density 残差 `−0.36`；过程中能量在约 `−297231` 至 `−1199` eV 间跳变，未进入稳定区。这些 partial 数值不作为标签。
- **诊断 1**：maxiter=240、Mixer β=0.10/nmaxold=8/weight=100；第 13 步人工停止，能量 `−145095.064419 eV`、残差 `+6.27/+0.12`。停止原因是仍有强振荡；没有到 maxiter，也没有墙钟或资源终止记录。
- **诊断 2**：maxiter=300、β=0.02/nmaxold=8/weight=100；第 34 步人工停止。能量曾短暂到 `−5.514891 eV`，随后退回 `−132.587053 eV`；末残差 `−1.55/−1.08`，仍远离完整收敛阈值，轨迹不是单调下降。不是到迭代上限或资源限制。
- **诊断 3**：maxiter=300、β=0.05/nmaxold=3/weight=100；第 123 步人工停止，末能量 `−134983.495505 eV`、残差 `+5.05/+0.04`，中间仍大幅跳变。不是到迭代上限或资源限制。

三次诊断中，β=0.10 直接放大早期振荡；β=0.02 一度降低残差但没有稳定能量；β=0.05 缩短历史到 3 项也未解决振荡。三次计算都使用相同 GPAW 默认 PPCG（日志为 niter=2），并从冻结结构启动；没有载入 checkpoint。三次实际 SCF 均由观察到不收敛趋势后人工停止，而非达到上限、墙钟超时、内存/CPU 限制或计算异常。启动阶段另有一次 `ModuleNotFoundError: No module named 'input'`，发生在首个 SCF 步之前；随后修正了副本 runner 的模块搜索路径，未把该启动错误计作 SCF 试算。

收敛日志均显示先“Diagonalizing LCAO Hamiltonian”，再“Converting LCAO to pw mode”；输入没有指定外部初始密度或 restart 文件。R1 收敛后记录了费米能级和占据表；R2 的未收敛日志没有最终 Fermi level/occupation summary，故不能判断 R2 是否发生了特定能带占据切换。没有发现 `*.gpw` 文件或可恢复的 R2 checkpoint。

## 可复用证据与数据边界

- 原 Phase2 commit：`0dd8628e4a91e88f31aba3fd7bbce1c8a40a2aa2`。原始 `phase2_8x8x1_20261012_exec01/R2/` 未修改。
- R1 仍为原有 SUCCESS：152 步，Energy `−188.92149225542852 eV`，force-consistent free energy `−189.16802546281826 eV`，forces/stress 已保存；本次未重算 R1。
- 原 R2 仍为 FAILED_SCF：160 步；Energy/free energy/forces/stress 均不可用。
- PAW 版本 `gpaw-data 1.2.1`；Ag、C、Si、Ti SHA 与固定参数锁一致。
- 三个诊断输出、状态 JSON、参数 runner、启动输出/错误回溯、原始副本结构和各级 SHA 清单均留在本目录；运行结果标为诊断，不进入 TRAIN/DEV/TEST。

## 唯一建议的下一判别点

请求 C 仅批准一个独立 R2 诊断 single point，固定 R2 坐标、cell、PBC、PBE/PAW、8×8×1、500 eV、Fermi–Dirac 0.10 eV、三个 1e-5 收敛阈值、原 Mixer β=0.05/nmaxold=8/weight=100、对称性、Poisson 参数、4 MPI ranks 和 OMP=1。

唯一电子求解器改动：把隐式默认 `PPCG(niter=2)` 显式设为 `PPCG(niter=4)`。已只读确认 GPAW 26.7.0 接受这一参数形式；未运行计算。理由是 R2 三次 mixer 试验都保留 PPCG niter=2，而原始 R2 的本征态末四步本征态残差在 +1.23 至 +2.12 间摆动、末值 +1.80，表明本征态更新不足与 density feedback 值得优先区分。代价是每个 SCF 步可能更慢、临时内存/计算开销上升；若真实原因主要在密度混合，该方案仍可能失败。不得据此自动启动第二点。

为避免将单纯增加上限误当修复，将 maxiter 从 160 限定提高到 200 作为唯一点的有限迭代预算；墙钟硬上限 4 小时。任何上限到达都标 `FAILED_SCF`，不重跑。runner 应为该副本设置每 40 步覆盖式保存一个诊断 GPW checkpoint（含 wavefunctions），保存完整 stdout/stderr 和 SCF log；启动前检查目标目录空间并校验 checkpoint 可读性。若无法保证有限墙钟停止后能写出有效 checkpoint，则保持 HOLD，不启动。

预期只判断“提高 PPCG 内部迭代是否稳定降低 R2 eigenstate/density 残差并完成原阈值”；这是一个诊断判别点，不保证恢复成功。收敛后才允许收集能量、free energy、forces、stress 并与 R1 的同配置结果比较。失败则停止并报告，不自动改 mixer、smearing、k 点、cutoff、结构或初始猜测。

执行状态：`EXECUTION_AUTHORIZED=false`。请 C 审阅并签发绑定本 proposal SHA、R2 input/structure SHA、参数锁 SHA 和 runner/collector SHA 的单点授权。当前保持 HOLD。
