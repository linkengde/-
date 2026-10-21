# Ag–Ti₃SiC₂ 三维双连续几何预览（V2）

这份结构用于查看用户参考图所示的三维互穿形貌。`USER_REFERENCE_STYLE_YZ_SLICE.png` 是单个 YZ 截面的相区示意：蓝色表示 Ag，灰色表示 Ti₃SiC₂；`BICONNECTIVE_3D_ATOM_VIEW.png` 是原子级斜视图。两张图展示同一类随机场形貌，不是热模拟轨迹。

## 结构信息

- 原子数：37,711；Ag 11,524，Ti 13,400，Si 4,112，C 8,675。
- Ag 质量分数：59.08 wt%。
- 晶胞：110 × 98.085 × 100.921 Å；边界为 X 非周期、Y/Z 周期。
- 材料的 X 坐标范围为 50.63–107.89 Å。
- 两相从独立的理想 FCC Ag 和 Ti₃SiC₂ 母晶格生成，再用固定随机种子的三维平滑随机场划分相区；没有复制历史熔池坐标。
- 界面几何去重后的最小异相距离为 2.200 Å。此项只表示几何距离，不表示势函数稳定。
- 连通性审计使用 Ag 3.15 Å、Ti₃SiC₂ 2.75 Å 的近邻阈值。两相最大连通分量均超过各自原子数的 99.8%，并在 Y/Z 周期方向贯穿；Ti₃SiC₂ 对近邻阈值敏感（2.55 Å 时主分量比例明显降低），因此该结果是几何拓扑检查，不是机械强度或动力学稳定性证明。

## 使用边界

**仅几何草模，未通过静态受力或动力学验证，不可把它当作已验证可加热模型。** 当前只用于浏览和讨论双连续相形貌；本目录没有启动 LAMMPS、TTM、DFT 或势函数训练。该结构也不是原 39,486 原子模型或已经运行的 766 原子加热试算，不能替换它们。

本小模型继承原参考结构的坐标偏移：材料下侧有较大 X 空区，而上侧仅留约 2.11 Å 到非周期盒边界。若后续要让原子从表面迁移/逸出或应用分区热边界，必须先重新布置 X 方向真空与支撑分区并重新审计；这不影响本次只看相形貌。

## 文件

- `RANDOM_BICONTINUOUS_V2.data`：LAMMPS atomic 格式几何。
- `RANDOM_BICONTINUOUS_V2.extxyz`：OVITO/ASE 可读结构。
- `PHASE_ATOM_GROUPS_V2.csv`：atom ID、元素、相、区域及随机场信息。
- `RANDOM_BICONTINUOUS_V2_AUDIT.json`：组成、晶格应变、距离和连通性审计。
- `INTERFACE_SITE_REMOVALS.csv`：界面重叠处理时被移除的初始位点。
- `build_bicontinuous_v2.py`、`build_random_bicontinuous.py`：几何生成脚本。
- `USER_REFERENCE_STYLE_YZ_SLICE.png`：单个相区切片预览。
- `BICONNECTIVE_3D_ATOM_VIEW.png`、`BICONNECTIVE_YZ_SLICES.png`：三维原子视图和多个切面。

完整文件校验值见 `SHA256SUMS_V2.txt`。
