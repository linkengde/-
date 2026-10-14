# 数据质量与独立性审查

## E2 根因审计状态

仓库中没有单独标名 `Stage E2` 的审计文件。此前 `potential_root_cause_final_assessment.md` 和 `TSCBOP_force_error_environment_analysis.md` 已完成组件力分解、实现/映射核验和空间误差分组；关键证据是 48 原子体相控制的 Ti/C 力误差已经很大，而 114 原子界面派生 TSC-only 对照总力 RMSE 为 2.5624 eV/Å。此前没有 TSCBOP 参数敏感性结果。本审计用隔离参数副本补了参数族方向敏感性，但不把它解释成唯一物理根因或正式拟合。

## DFT 方法与完整性

48 原子数据均为 GPAW 26.7.0 / ASE 3.29.0 / gpaw-data 1.2.1，PBE/PAW、500 eV 平面波、Fermi–Dirac 0.10 eV、SCF energy/density/eigenstates `1e-5`、4 MPI rank、PBC TTT。基线及 Ti/C 位移计算输出 `PW_PBE_forces`，单位 eV/Å；native extrapolated total energy 与 force-consistent free energy分别保存，单位 eV/cell。11 个结果档案、纯相包清单均通过 SHA 校验，SCF verification 为 PASS。

48 原子基线是 CIF-derived 固定晶胞近父结构，组成 Ti24Si8C16，cell 为 `[[6.115,0,0],[-3.0575,5.295745344,0],[0,0,17.6235]] Å`，并非已弛豫到零力的 PBE 平衡晶体：10×10×4 下 Fmax=0.08786 eV/Å、RMS=0.05699 eV/Å。原始 CIF 的可审计出处没有随该输入包完整归档。该点是已收敛固定几何参考，不应命名为“平衡体相”。

7 个基线 k 网格输出的 `input_geometry_sha256` 完全相同 (`157c1dab…7805c3e4`)；它们是数值采样对照，不是 7 个独立结构。参数敏感性只选用 10×10×4 作为基线，避免重复 K 网格加权。四个响应输入均为同一 48 原子父结构上的单原子 z 位移 ±0.02 Å，Ti 4f local ID 3、C 4f local ID 9，使用相同 10×10×4 方法；它们是一个 Ti 位移族和一个 C 位移族，不是独立来源。DFT 对称力响应分别为 Ti 17.4712、C 15.2505 eV/Å²。没有 Si 位移 DFT。

48 原子基线及四个位移的 extxyz、summary 和 gpaw.log 中均没有应力张量记录；因此这批数据不能约束体弹性或剪切响应。native/free energies可用于一致约定下的相对能量目标，但不能与未校零的 LAMMPS 绝对能量直接相减。

## 表面/界面记录

有一份 114 原子 `Ti57Si19C38` TSC-only 固定几何 DFT：从已优化的 135 原子 Ag21Ti57Si19C38 界面移除 Ag，保留其余坐标、cell 和 PBC。子集输入 SHA=`260ae9b736cb95a641b4d1eae9ffbec5711f2ffe895ed1f0acc7e25ff637fdf1`，父优化结构 SHA=`eda2586c63c6c3c793ffce99417e7c04569941b0a2a916fda1b5402b7ba52d9f`。PBE/PAW 500 eV、0.10 eV smearing、Γ 1×1×1、1e-5，SCF 从有效 partial.gpw 恢复并在44步收敛；能量 -872.877543 eV/cell、free energy -874.741340 eV/cell，114条力 eV/Å，应力 eV/Å³可读。组件清单逐项 SHA 校验通过。

该 slab 不是独立构造并弛豫的 clean Ti₃SiC₂ 表面；它是移除 Ag 后的界面派生冻结 slab，最大骨架力约0.9757 eV/Å，不能作表面平衡能或唯一表面模型。只有这一份 slab 家族，不构成表面独立验证。另有 135 原子全界面优化结果（Fmax=0.04278 eV/Å），但含 Ag；它是交叉界面背景，不能作为 TSC-only 骨架标签替代。

## 角色冻结与独立性

所有原记录原封不动保留其已有角色：48 原子及位移仍为 `numerical_pure_phase_control_not_training_or_independent_validation`；114 原子 slab 仍为单个分项诊断；135 原子全界面只作上下文。`role_freeze.json` 将当前 TRAIN、DEV、TEST 冻结为空，防止把控制/数值对照擅自并入数据集。若经后续审核将一个 parent family指定为探索性 TRAIN，必须另建 split manifest；当前没有独立 DEV 或 TEST，任何拟合只能被称为探索性，不能声称泛化。

结论：当前有一个48原子体相父族（基线+Ti/C响应）和一个界面派生slab族；K网格对照不增加几何独立性。没有已验证平衡体相应力/弹性响应、Si局部扰动或独立clean-surface验证。
