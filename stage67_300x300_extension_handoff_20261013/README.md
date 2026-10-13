# Stage67 300×300 Å 扩展建模交接包

状态：HANDOFF_REFERENCE_PACKAGE；不是可直接运行的 300×300 Å 模拟输入包。

## 目标与来源

- 目标：以 Stage67 已审计的 Ag–Ti₃SiC₂ 熔池作为结构来源，准备扩展 Y/Z 横向尺寸至约 300×300 Å，后续用于压力/熔坑研究。
- 来源帧：Stage67 200 Å 横向模型，step 1650，t = 0.165 ps，303,372 atoms，盒子 150.1×200×200 Å，边界 p p p。
- 随包包含原始 step1650 dump、对应 200 Å checkpoint.restart 和 checkpoint_CellTl_state.dat，用于来源核验与后续几何/场映射设计。
- 已有稳定熔池报告将 0.155–0.165 ps（step1550–1650）定义为采样稳定窗口。证据范围仅为该 10 fs 窗口，不表示长期/热力学稳定。

## 重要限制：目前不能直接跑 300×300 Å MD

本包没有已经建立并验证的 Stage67 300×300 扩展结构，也没有与该目标几何匹配的可运行 LAMMPS 输入和 15×30×30 TTM 空间映射表。因此，包名表示交接目标，不表示 300 Å 模型已生成。

扩展前必须由云端明确并记录：
1. 从 200×200 到 300×300 的确定性几何生成规则（200→300 不是整数倍平铺；不得简单缩放坐标或静默改变原子密度、Ag/Ti₃SiC₂比例、枝晶形貌或自由表面）。
2. 新结构原子数、成分比例、周期缝合/边界连续性及结构检查。
3. 若保持约 10 Å 网格，按目标尺寸重建 15×30×30 的 cell_material_fraction、G_mix、初始 Te 索引和热源分布。input/ttm/200A_spatial_maps_NOT_FOR_300A/ 仅供来源追溯，禁止直接作为新网格输入。
4. 与新结构一致的 LAMMPS 输入/脚本及其引用路径。Python driver 是参考材料且含本机路径；Stage67 manifest 引用的 in.regions_x_surface.inc 当前未在来源目录找到，不能用其他模型的 region include 替代。
5. 在任何正式长跑前，核实势函数、TTM实现/单位和热源映射；交叉势 potential.interface_calibration_v0.inc 明确是 calibration v0，未获生产级验证。

## 保持项

来源资料记录旧 Stage67 设置为：Ag EAM、Ti–Si–C Tersoff 加 calibration-v0 交叉项；MD timestep 0.1 fs；约 10 Å TTM cell；周期边界类型 p p p；热源使用旧 Stage67 的空间审计表。新几何中要保持物理定义和边界类型，须在重新映射后核验；旧 CSV 网格索引不可原样照搬。

液态判据引用现有 step1650 方法：Ag type 4、PTM Structure Type=Other、CellTl_valid=1 且 CellTl≥975 K；最大簇 cutoff=3.74625 Å。JSON 还记录 Q6 threshold 0.446766，但现有包没有所有 Ag 的逐原子 PTM/Q6标签或完整 LiquidCandidate列表，不能据汇总文件重新生成全体液态原子集合。

## 文件说明

- input/source_200A_step1650/：精确来源帧、restart、TTM侧车及其温度/能量审计。侧车是旧 200 Å网格状态，不可直接赋给 300 Å网格。
- input/potentials/：Stage67旧相内势和交叉势文件；交叉势是未验证的 calibration-v0。
- input/ttm/material_tables/：材料性质表，部分温度范围含诊断/外推；核对单位并确认实现匹配后再用。
- input/ttm/200A_spatial_maps_NOT_FOR_300A/：旧网格的材料分数、耦合、初始Te和热源映射。
- input/implementation_reference/、input/analysis_scripts/：带路径/版本限制的参考脚本，不是独立可运行流程。
- reports/：稳定窗口、PBC、TTM、势函数、热源和既有 300 Å准备度审查。
- INVENTORY.csv：每个已打包来源文件的路径、大小、SHA256和适用状态。
- reports/SHA256SUMS.txt：对本包其余文件的逐文件SHA256（不包含本清单自身）。

## 源状态的冲突记录

原始 checkpoint 的 checkpoint_cluster_summary.csv/checkpoint_status.md 记录 liquid-like Ag 数为 0，而独立稳定熔池报告和 step1650 Q6 JSON 记录 8,728 个候选 Ag、最大簇 1,032 个。前者以原始文件保留用于追溯，但不作为熔池判定；云端应先解释该差异，不能把两套统计合并或静默覆盖。

## 操作边界与建议

准备本包仅复制原始文件并新增说明/清单；未修改源文件、未运行 LAMMPS/MPI、未推进 MD。云端第一步建议只读验证文件哈希、盒尺寸/原子数、restart/TTM侧车兼容性和待补项；完成几何规则与新网格映射设计前不要启动长时间模拟。
