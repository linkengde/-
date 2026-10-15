# LOW_MISMATCH_135_DFT_LITE_SINGLE_POINT_REPORT

## 结构与参数核验

- 结构文件：`Ag_Ti3SiC2_12p5198deg_135atom_dft_lite.extxyz`
- SHA-256：`1a3e341e1132bba9b349ab0fe6e3a3be11d1e1661671b96f6810a54f2e1f54ae`（与任务要求一致）
- 来源：fresh-fetch 后由提交 `f1a080e14cd5ffd9739dcd2a3df2d188c853a8a8` 中的 217 原子结构按既有完整层裁剪流程确定性生成；217 原子源 SHA-256：`945c5e214c1ceaeba4125ff785f4092e14858c75a985061a7841da3e5ce6a9a2`
- 原子数/组成：135；Ag21Ti57Si19C38
- PBC：`[True, True, True]`；晶胞（Å）：`[[13.38617865561341, 0.0, 0.0], [6.693089327806705, 11.592770775358236, 0.0], [0.0, 0.0, 37.9451596]]`
- z 真空：29.10793160 Å
- 扭转角：12.519829797°；Ag 面内应变：1.102898177%；Ti₃SiC₂ 面内应变：0.0%
- PAW Ag/Ti/Si/C SHA-256 已核验并绑定在 `single_point_gamma_1x1x1/INPUT_PREFLIGHT.json`。
- GPAW 26.7.0、ASE 3.29.0、OpenMPI 5.0.7；4 MPI ranks、OMP=1。
- 科学参数：PBE/PAW，500 eV，Fermi–Dirac 0.10 eV；固定几何；energy/density/eigenstates 收敛阈值均为 1e-5。

## 1×1×1 gamma-centered 固定几何单点

- 状态：**CONVERGED**，44 SCF steps。
- 总能：**-932.134316825458 eV**。
- 力一致自由能：**-933.557075245305 eV**。
- 力：135 个原子均有力；最大力 **7.216193 eV/Å**。
- 应力（eV/Å³）：`[[0.0004006463448134574, 4.874056252199953e-07, 9.653586257309315e-07], [4.874056252199953e-07, 0.00039959681914613656, 5.751014263512344e-08], [9.653586257309315e-07, 5.751014263512344e-08, -0.013150069234390858]]`。
- 界面力定义：最高 z 的 Si 层与最低 z 的 Ag 层，共 40 个原子；界面最大力 **7.216193 eV/Å**，RMS **2.968046 eV/Å**。
- 最近 Ag–Si 对：Ag index 114 与 Si index 3，距离 2.000 Å；两者力模分别 7.216193 和 7.014359 eV/Å，主要沿 z 方向相反。
- 峰值 cgroup 内存：13,969,739,776 bytes = 13.01 GiB；基线 5.32 GiB；墙钟 3793.5 s。
- 计算过程中坐标、晶胞和 PBC 保持不变；无结构优化。

固定几何下界面力较大，特别是 2.0 Å 最近 Ag–Si 对的相反 z 向力。该结果只证明本配置完成了 SCF 单点；不代表此结构达到力学平衡，也不授权后续优化。

## 2×2×1 评估

- 参数保持一致，仅将 gamma-centered k 网格改为 2×2×1。
- 状态：**ABORTED_MEMORY_GUARD_BEFORE_SCF**；GPAW 初始化期间达到 28 GiB 安全线，SCF 尚未开始。
- 峰值采样：30,111,612,928 bytes = 28.04 GiB；cgroup 上限 32.00 GiB；安全线 28.00 GiB。
- 作业在 762.7 s 后被内存保护器停止。未产生可用能量、力或应力；无自动重试、无参数更改。
- 结论：当前 32 GiB cgroup 环境下，2×2×1 未能通过内存安全预检，应保持 HOLD；1×1×1 结果有效且已完成。

## 归档文件

- 1×1 结果、SCF 日志、全原子力和界面力：`single_point_gamma_1x1x1/`
- 2×2 停止证据：`single_point_gamma_2x2x1/RESOURCE_USAGE.json`、`FAILURE.json`、`gpaw.log`
- SHA 清单：`SHA256SUMS.txt`
