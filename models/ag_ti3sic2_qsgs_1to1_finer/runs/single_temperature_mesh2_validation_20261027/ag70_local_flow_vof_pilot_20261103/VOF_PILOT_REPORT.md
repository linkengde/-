# 液态流动求解链路试验：Poiseuille 通道与 VOF 静滴

## 目的和范围

为“液态 Ag 流动—自由表面变形—熔坑/喷溅”研究先验证独立流动求解链路。使用 OpenFOAM 2406 的 `simpleFoam` 对单相直通道压差流动与解析 Poiseuille 解比较，再用 `interFoam` 做无重力、无外部载荷的二维圆柱液滴静态表面张力试验。全部均为**无量纲数值方法验证**，不是 Ag 物性算例，不是电弧计算，也不代表熔池或喷溅结果；没有使用 MD。

## 压差通道基准

120×40×1 单相层流网格，长宽比 12:1；运动黏度 `nu=0.01`、无量纲压差 `Δp=1`。该压差仅作为解析解验证驱动，不是实际电弧压力。对收敛后的中心段速度剖面与 `u(y)=Δp·y(H−y)/(2 nu L)` 比较：

| 指标 | 结果 |
|---|---:|
| 解析最大速度 | 1.0416667 |
| 模拟最大速度相对误差 | −9.3×10⁻⁹ |
| 速度剖面归一化 L2 误差 | 0.0625% |
| 单位跨度流量误差 | 0.125% |
| 最后写出间隔剖面变化 | 4.8×10⁻⁹ |
| 累计连续性误差 | −5.61×10⁻⁷ |

这个结果支持压力—动量离散能正确复现单相通道流；它不验证气液界面、润湿、骨架阻碍或 Ag 物性。

前两组保持相同域、液滴半径和无量纲参数，只加密网格；第三组单独提高无量纲黏度并延长计算时间，用于检查瞬态阻尼：

- 80×80×1，`dx=0.025`；
- 160×160×1，`dx=0.0125`；
- 液相、气相均设 `rho=1`、`nu=0.01`，`sigma=1`，`g=0`，液滴名义半径 `R=0.4`；
- 网格单元中心二值初始化液相，计算至无量纲 `t≈0.2`。
- 阻尼检查使用 `nu=1`，其余参数相同，计算到 `t≈0.5`。

二维圆柱界面的解析压差为 `Δp=σ/R_eq`，其中 `R_eq` 按数值守恒的液相面积求得。

## 结果

| 网格 / `nu` / 终止时间 | 面积相对漂移 | 解析压差 | 数值液/气体相平均压差 | 压差误差 | 最终瞬时最大速度 | 轮廓长短轴比 |
|---|---:|---:|---:|---:|---:|---:|
| 80×80 / 0.01 / 0.2 | `5.03×10⁻¹³` | 2.4880 | 2.1075 | −15.30% | 0.2816 | 1.0000 |
| 160×160 / 0.01 / 0.2 | `−1.51×10⁻¹¹` | 2.4957 | 2.1674 | −13.16% | 0.2979 | 1.0000 |
| 80×80 / 1 / 0.5（阻尼检查） | `−6.80×10⁻¹¹` | 2.4880 | 1.9363 | −22.18% | 0.0111 | 1.0000 |

累计连续性误差约 `10⁻²⁰`，液相面积守恒良好，液滴质心居中且轮廓保持圆形。前两组加密一倍后压差误差只从 15.3% 降到 13.2%，没有显示充分收敛；末时刻速度仍约为 0.28–0.30。阻尼检查把速度降到 0.011，但压差误差升至 22.2%，压力跳跃从 `t≈0.02` 的 1.816 持续升至 `t≈0.50` 的 1.936，仍未收敛到解析值。不能单凭“速度变小”判定达到静态平衡。结论是：**单相通道动量基准通过；VOF 求解器能运行且液相质量守恒良好，但表面张力/曲率平衡验证未通过，不能进入 Ag 熔池自由表面预测。**需要先调查 VOF 曲率离散、压力跳跃评估和网格/时间收敛，再重复静滴及毛细波/润湿测试。

各保存时刻的面积、压力跳跃和最大速度见 `static_drop_time_history_80.csv`、`static_drop_time_history_160.csv`、`static_drop_time_history_viscous80.csv`；汇总见 `mesh_sensitivity.csv`。

## 与 Ag 熔池模型的关系

本试验不含 Ag 真实密度、黏度、表面张力，不含 Ti₃SiC₂ 骨架、焓法、传热、界面润湿、蒸发、反冲压力或电弧载荷。现有 64×64×50 焓法数据只提供温度、液相分数和骨架掩膜；它不能直接提供自由表面流动条件。即使静滴基准通过，后续也要先建立更大且包含气相的物理域，并用实验/文献校准液态 Ag 与电弧边界参数。熔坑喷溅阶段还必须加入自由表面和气相耦合；若要预测蒸发喷溅，需要再加入质量损失与蒸气反冲模型。

## 复现

使用的容器镜像为 `opencfd/openfoam-default:2406`，下载时解析出的镜像 digest 为 `sha256:dd5aa20630a55722663bf83ba0cb74870cba130081303e32e3865007fa2aa35a`。运行环境需安装 Docker 并已有此镜像。

```sh
python prepare_poiseuille.py
./run_poiseuille.sh
python analyze_poiseuille.py

python prepare_case.py --resolution 80
./run_case.sh static_drop_case
python analyze_case.py --resolution 80

python prepare_case.py --resolution 160
./run_case.sh static_drop_case_160
python analyze_case.py --resolution 160

python prepare_case.py --resolution 80 --case-name static_drop_case_viscous80 --nu 1 --end-time 0.5 --max-dt 0.0005
./run_case.sh static_drop_case_viscous80
python analyze_case.py --resolution 80 --case-name static_drop_case_viscous80 --output-tag viscous80
```

`poiseuille_metrics.json`、`poiseuille_profile.csv` 和 `poiseuille_profile.png` 保存解析解对照；`static_drop_metrics_*.json` 是 VOF 指标；`static_drop_benchmark_*.png` 显示液相轮廓和末时刻速度；`static_drop_final_fields_*.npz` 保存末时刻液相分数、压力、速度和网格坐标。完整输入字典保存在 `poiseuille_case` 和三个 `static_drop_case*` 目录中；CSV 保留各写出时刻的统计，不需要把全部中间场上传到仓库。
