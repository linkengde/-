# QSGS 更细骨架候选 05

状态：`UNRELAXED_GEOMETRY_CANDIDATE_NOT_FOR_HEATING`
日期：2026-11-25

## 建模方法

本候选沿用原随机 QSGS 生长流程，没有手工绘制或摆放骨架枝条。网格为 `100 × 100 × 75`、体素尺寸 `2 Å`，随机种子 `20261125`，核心成核概率 `p₀=0.064`，26 邻域生长概率为 `0.18/r`。自动后处理按既有流程处理体素连通性、原子最近邻间距和两相质量配比。

为让本次较细的骨架在 Ag 连通性筛选后仍有足够 Ag 原子达到 1:1 质量比，质量平衡搜索的预筛余量阈值由候选 03 的 `1.08` 调为 `1.20`。用旧阈值试算时，Ag 连通筛选后为 62,926 个原子，低于所需的 66,686 个，构建失败；采用 `1.20` 后构建成功。该调整只影响自动形貌筛选余量，没有手工改变枝条或局部坐标。

## 与候选 03 的形貌比较

两候选使用相同的 `200 × 200 × 150 Å` 材料区、`200 × 200 × 300 Å` 模拟盒、x/y 周期边界、Ti₃SiC₂ 晶格取向、Ag 晶格常数和约 1:1 两相质量比。候选 05 将 QSGS 核心成核概率由 0.024 提至 0.064，并使用新随机种子。

| 指标 | 候选 03 | 候选 05 |
|---|---:|---:|
| QSGS 成核概率 | 0.024 | 0.064 |
| 骨架体素比例 | 0.472907 | 0.425216 |
| 三轴平均骨架线截长度 | 5.752 Å | 4.513 Å |
| 相对候选 03 的线截变化 | — | 减少 21.5% |
| Ti₃SiC₂ 原子数 | 236,880 | 212,424 |
| Ag 原子数 | 71,630 | 64,234 |
| 总原子数 | 308,510 | 276,658 |
| Ag:Ti₃SiC₂ 质量比 | 1.000006:1 | 0.999994:1 |

线截长度是体素网格上沿 x、y、z 三轴统计的平均连续骨架长度，只用于比较统计形貌，不能解释为每条原子级连接颈的厚度或机械强度。QSGS 生成的是随机双连续多孔结构，不保证形成规则树枝状外观。

## 几何检查

- Ti:Si:C 原子数严格为 3:1:2。
- Ag 和 Ti₃SiC₂ 在原子近邻图中各为一个连通分量，并分别跨越 x/y 周期接缝和材料 z 厚度。
- 两相质量比为 `0.9999942:1`。
- 最近的 Ag–Ti、Ag–Si、Ag–C 间距分别约为 2.228、2.228、2.228 Å；这是建模几何筛查，不代表势函数稳定性或界面能验证。
- 材料区为 `z=0–150 Å`，上方真空 150 Å；底部 `z=0–20 Å` 标记为固定层。

## OVITO 图

- `thinner_scaffold_OVITO.png`：完整模型和中心剖切双视图。
- `thinner_scaffold_OVITO_full.png`：完整模型视图。
- `thinner_scaffold_OVITO_center_cut.png`：中心剖切视图。

蓝色表示 Ti₃SiC₂ 骨架，浅灰色表示 Ag。图片由 OVITO 渲染脚本输出。

## 复现

在仓库根目录、且依赖已安装时运行：

```bash
python models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/geometry_constraint_B/thinner_scaffold_candidate_05_20261125/build_thinner_scaffold.py
python models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/geometry_constraint_B/thinner_scaffold_candidate_05_20261125/render_thinner_scaffold_OVITO.py
```

构建依赖列于 `models/ag_ti3sic2_qsgs_1to1_finer/requirements-build.txt`，OVITO 可视化依赖列于 `models/ag_ti3sic2_qsgs_1to1_finer/requirements-visualization.txt`。`geometry_report.json` 记录本次完整参数、计数、晶格来源、间距和连通性统计；`SHA256SUMS.txt` 用于校验文件。

## 使用限制

这是新的未弛豫几何候选，不是之前 25.5 ps 轨迹的续算结构，也不能复用旧 restart。当前没有运行势能最小化、300 K 冷态稳定性测试、MD、加热或相变验证。继续热模拟前应先用选定势函数做静态受力和短时冷态检查。候选 05 的统计线截较候选 03 更小，但仍是 QSGS 随机多孔骨架；是否满足最终视觉和力学需求需在 OVITO 视图及后续验证中判断。
