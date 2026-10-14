# 最小补充 DFT 计划（仅规划，未执行）

## 现有数据清点

以下结果可以复用，不重复计算：

1. **48 原子体相。** Ti₃SiC₂ PBE/PAW 固定结构存在多个 k 网格结果（包括 4×4×2 至 10×10×4）。这些结果主要检查数值采样，父几何相同，不能算成多个独立构型。TSCBOP 对其中 k4/k6 同结构 force RMSE 为 1.9965/1.9872 eV/Å。
2. **局部 Ti/C 位移响应。** 已完成 4 个固定晶胞 DFT 点：Ti 4f 与 C 4f 各自沿 z ±0.02 Å，10×10×4。它们共享一个父结构，不是独立验证集。
3. **界面派生的 TSC-only slab。** 已有 114 原子 Ti₅₇Si₁₉C₃₈ 冻结坐标单点 DFT，含底部 Ti 表面及上部 Si 终止，来自去除 Ag 的已弛豫界面构型。可复用为表面/大力诊断，但不代表独立松弛干净表面或表面能。
4. **完整 Ag 界面。** 已有一个 135 原子 registry 的 DFT 优化路径和最终构型。该单一轨迹父族不足以支持跨 registry cross 拟合。
5. **未执行提案。** Ti₃SiC₂ 体积 ±2% 和剪切 ±0.5% 的四个 `PROPOSAL.extxyz` 已存在，但没有相应 DFT 结果；不把提案当作标签。

纯相候选体积/剪切提案目录：`research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/inputs/`。DFT 位移执行记录：`research/mace-v12-transfer/coordination/reports/window-b/v19_Ti3SiC2_backbone_force_response_execution_window_b/FINAL_RESULTS.md`。本节只依据归档文件清点，没有打开封存数据或改变数据角色。

## 优先补充：先判定骨架项，不拟合

1. **先复用已有 baseline 和 4 个位移点。** 不再重算 k 网格或 ±0.02 Å 位移。
2. **下一批最小 bulk response 为四个固定几何单点：** 体积 −2%、+2%；面内剪切 −0.5%、+0.5%。用同一 48 原子父结构、PBE/PAW、10×10×4、相同 cutoff/展宽/收敛标准，输出 total/free energy、forces、stress；原子与晶胞变化按 proposal manifest 冻结。它们直接检验 TSCBOP 体积/剪切响应和 DFT 曲率方向。
3. **资源预算：** 已完成的相同 48 原子基准/扰动 DFT 点约 2.83–4.20 小时/点（4 MPI ranks）。按此经验，四点顺序执行约 **11–17 小时墙钟**；这是历史时间范围外推，不是本轮新测时。当前环境可做单点，但仍需后续正式执行授权。
4. **表面门：** 先将现有 114 原子 TSC-only DFT 力标签纳入诊断，不另造重复计算。若拟合/应用明确要求正确的干净表面弛豫或表面能，再由 C 冻结一份有可靠结构来源的 6 层、Ti/Si 两终止 slab 几何；至少做一份固定几何 DFT force 单点及一个确定性表面扰动点。现有 interface-derived slab 不能替代这项独立检验；没有结构来源审定前不启动。

这六类信息（bulk baseline、四个现有位移、四个新体积/剪切、现有 slab）足以做“能否继续保留 TSCBOP 形式”的小规模判断，不足以证明温度、缺陷和所有终止面可靠，更不等于已拟合。

## Ag–Ti/Si/C cross 的最小后续数据

保留当前已归档有效力分解，不改 TSCBOP 或 Morse。本轮已有单个 135 原子 interface geometry，只能作基线。若骨架模型先通过 bulk/surface gate，后续最小 cross screening 建议由 C 冻结：

- 在 Si 终止面上增加 **两个非等价横向注册**；
- 每个注册三个法向间距（压缩、近平衡、拉开），间距覆盖当前最低附近约 ±0.15–0.30 Å，形成 **6 个 full-interface fixed-geometry DFT 单点**；
- 再增加 **2 个确定性小位移**点测试近平衡力响应；
- 加 **1 个大于 8 Å cross cutoff 且无周期镜像相互作用**的分离参考。

这是约 9 个新增 full-system 构型的起始筛选，不是最终可迁移数据集。full DFT energy/forces/stress 是统一 MACE 训练的标签；若仍要把 cross 力从骨架中分离，则在少量代表 gap/registry 上补同几何 Ag-only 与 Ti₃SiC₂-only 参考。不得用 cross Morse 去抵消 TSCBOP 骨架力误差。所有 TRAIN/DEV/TEST 角色由 C 在执行前冻结，当前没有结构生成、DFT 或拟合授权。
