# Stage67 300×300 Å 扩展交接前置审查

## 判定

PACKAGE_READY_FOR_300A_MD = NO
PACKAGE_READY_FOR_CLOUD_REVIEW_AND_GEOMETRY_DESIGN = YES

## 已提供且可追溯

- Stage67 step1650 200 Å横向原子帧：303,372 atoms，150.1×200×200 Å，p p p。
- 同时刻来源 checkpoint restart、CellTl/Te/U网格状态侧车及温度/能量历史。
- Stage67 Ag EAM、Ti-Si-C Tersoff、Ag-Ti/Si/C calibration-v0 include。
- Ag/Ti₃SiC₂ Ce、kappa、G材料表、旧200 Å空间映射/热源审计。
- step1650 Q6/PTM汇总、最大簇ID和0.155–0.165 ps稳定窗口报告。
- 旧driver、分析定义以及PBC/尺寸审计参考。

## 仍缺少/未通过

1. 300×300扩展坐标文件：不存在可证明从Stage67 step1650扩展而来的300×300 data或restart。
2. 扩展生成规则：200→300不是整数倍周期平铺；缺少经确认的枝晶延拓/边界缝合规则。不能把Stage65 300 Å冷态候选冒充为该模型。
3. 300 Å TTM网格：旧表是15×20×20。若目标仍约10 Å网格，应重新映射到约15×30×30。新材料分数、G_mix、初态Te和热源空间表缺失。
4. 可直接运行的模型输入：没有与step1650稳定状态及目标300 Å几何匹配的完整LAMMPS input/启动脚本。现有Python driver和分析脚本是旧路径/部分实现的参考。
5. Stage67 inputs manifest引用的in.regions_x_surface.inc未在当前源路径找到；不以其他模型的region include替代。
6. 物理资格：交叉势calibration-v0未获生产级验证；部分Ti₃SiC₂ kappa/G/Ce是诊断或外推数据。正式长跑前需独立核验。
7. 液态逐原子数据：现有JSON是聚合结果，缺少全体候选Ag的逐原子PTM/Q6标签和完整LiquidCandidate列表。

## 关键不兼容性

| 项目 | 来源状态 | 300×300要求 | 处理 |
|---|---:|---:|---|
| 横向盒长 | 200×200 Å | 300×300 Å | 明确几何拓展规则后重建 |
| TTM网格 | 15×20×20 | 约15×30×30 | 重新映射，不复用旧空间索引 |
| restart | 200 Å step1650 | 新盒子/新原子集合 | 只作源状态，不得直接续跑 |
| CellTe/U/Tl sidecar | 200 Å旧网格 | 新网格 | 需重新映射并验证守恒 |
| 边界 | p p p | 保持周期边界类型 | 验证扩展边缘的结构连续性和密度 |

## 源数据差异

checkpoint status/cluster summary的候选液态计数为0，而稳定窗口报告与step1650 Q6 JSON提供8728候选Ag、最大簇1032。此差异已保留并明确标注。不要用checkpoint摘要的0覆盖稳定窗口判据，也不要将聚合JSON误认为完整逐原子标签。
