# Candidate 05：3.0 ps 局部加热阶段记录

日期：2026-11-26
状态：加热阶段完成；3.0 ps 检查点可续算；尚未开始冷却。

## 模型和运行条件

- QSGS candidate 05，未弛豫初始几何；盒子 200 × 200 × 300 Å，材料占 z=0–150 Å，x/y 周期、z 非周期。
- 共 276,658 个原子：Ag 64,234 个，Ti/Si/C 骨架 212,424 个；可动 Ag 54,540 个，初始 z<20 Å 的 Ag 固定。
- 可动 Ag 初始速度按 300 K Maxwell 分布生成。LAMMPS 默认 `thermo temp` 把固定骨架等原子也计入自由度，所以初始 run0 记录的 59.14 K 不是可动 Ag 温度；后续温度均用 TTM 耦合 Ag 原子或热源加权统计。
- Ag–Ag 使用 Williams–Mishin–Hamilton Ag EAM（`potentials/Ag_u3.eam`）。Ti/Si/C 全部固定，Ag–骨架使用元素相关的 WCA/LJ 排斥接触模型。
- 该骨架在本轮只提供几何支撑，不振动、不传导晶格热；Ag–骨架也不是经 DFT 验证的真实界面。此模型是探索性理想边界模型。
- 热源沿用先前的 XY 圆对称高斯分布（R=30 Å）乘 z 向表面指数衰减（skin=20 Å），峰值 5 eV/(ps·Å²)，时间步长 0.5 fs。使用仓库中的 Python TTM 参考耦合器。
- 从 300 K 速度初始化开始，热源连续运行到 3.0 ps。总计 6,000 步；每 0.5 ps 保存一帧及原子/TTM 检查点。

## 运行结果

| 量 | 3.0 ps 结果 |
|---|---:|
| TTM 耦合 Ag 原子平均晶格温度 | 1,387.23 K |
| 热源加权局部晶格温度 | 1,388.53 K |
| 电子网格最高温度 | 6,479.65 K |
| 累计热源输入 | 8,958.60 eV |
| 累计电子到晶格能量转移 | 1,581.21 eV |
| 能量账本残差 | −0.0808 eV |
| 可动 Ag 最大力（本段历史最大） | 16.6659 eV/Å |
| LAMMPS 安全停止阈值 | 50 eV/Å |

所有加热段都正常完成，LAMMPS 未报告原子丢失或超过 50 eV/Å 的安全停止条件。最大力仍偏高，结果只用于探索性观察。

## 局部结构和迁移证据

热区按中心 r_xy≤30 Å、z≥120 Å 定义，初始包含 1,043 个 Ag 原子。对同一批原子，0→3 ps：

| 指标 | 0 ps | 3 ps |
|---|---:|---:|
| q6 平均值 | 0.6673 | 0.4974 |
| qbar6 平均值 | 0.6097 | 0.2714 |
| qbar6>0.45 的比例（仅作相对筛查） | 100% | 4.31% |
| 去除热区质心漂移后的 MSD | 0 Å² | 11.175 Å² |
| 对应 RMS 位移 | 0 Å | 3.343 Å |

局部 q6/qbar6 明显下降，MSD 持续增加，且热区温度超过 Ag EAM 参考熔化温度范围。这些指标共同支持“热区 Ag 强烈失序并呈液态化行为”的探索性判断。它们还不足以证明达到热力学平衡熔化。

PTM 在初始帧中已把约 99.1% 热区 Ag 判为 `Other`，3 ps 时为 100%；因此这次 PTM 分类几乎没有区分能力，不能把 PTM `Other` 单独解释为液态。qbar6 的 0.45 切分也仅用于同一条轨迹的相对趋势，不是通用相界阈值。

## 图像与数据

- `candidate05_heated_3ps_OVITO_data_view.png`：OVITO 读取 3.0 ps dump 并计算 Ag-only qbar6，再由 Matplotlib 软件后端绘制的三维视图；红色表示较低 qbar6，蓝色骨架为固定支撑。
- `continued_2p5_to_3p0ps/frame_03.000ps_step06000.dump.gz`：3.0 ps 完整粒子坐标、速度和受力。
- `continued_2p5_to_3p0ps/lammps.restart` 与 `ttm_state.npz`：成对续算检查点；配套 `CHECKPOINT.json` 含校验信息。
- `analysis_q6_ptm_msd_0_to_3ps/q6_ptm_msd_timeseries.csv`：已更新到 3.0 ps 的全 Ag 与热区结构指标。
- `analysis_q6_ptm_msd_0_to_3ps/msd_timeseries.csv`：已更新到 3.0 ps 的 MSD 数据。

OVITO 的 Tachyon/标准渲染在本环境无显示会话下仍被 Qt offscreen Vulkan 集成阻挡；因此 PNG 是基于 OVITO 导出的粒子/qbar6 数组用软件后端绘制，不冒称为 OVITO 原生 Tachyon 渲染图。

## 判定与续算

- **加热运行：探索性通过。** 3 ps 检查点有效，热区表现出明显结构失序和原子迁移。
- **真实 Ag/Ti₃SiC₂ 触头预测：未验证。** 刚性骨架不参与晶格传热，Ag–骨架作用是人为设定的排斥边界。
- **凝固：未测试。** 热源仍处于开启段，当前结果不能说明停止加热后会凝固。
- **后续续算：** 使用 `continue_heat_candidate05.py`，将 `AG_TTM_CHECKPOINT_DIR` 指向 `continued_2p5_to_3p0ps`，将 `AG_TTM_CONTINUE_OUTPUT_DIR` 设为新的输出目录，并设置 `AG_TTM_CONTINUE_STEPS=1000` 可按相同参数续算 0.5 ps。续算前不要覆盖原检查点。

复现输入：`ttm_heat_trial_candidate05.py`、`continue_heat_candidate05.py`、`cold_stability_0p05ps.py`；势和模型输入的哈希由各检查点 `CHECKPOINT.json` 记录。此候选几何没有经过弛豫，300 K 运行前仅通过短时稳定性筛查。
