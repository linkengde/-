# 液态流动求解链路试验：通道流与 VOF 静滴

## 范围

为后续研究液态 Ag 流动、自由表面变形和熔坑/喷溅，先验证独立流动求解链路。使用 OpenFOAM 2406：`simpleFoam` 做单相压差通道流解析解对照；`interFoam` 与 `interIsoFoam` 做二维圆柱液滴表面张力平衡测试。全部是无量纲数值验证，不是 Ag 物性或电弧算例，没有使用 MD，也没有改动生产焓法求解器。

## 压差通道基准

120×40×1 单相层流网格，通道长宽比 12:1，运动黏度 `nu=0.01`，施加无量纲压差 `Δp=1`。该压差只用于解析解验证，不代表真实电弧压力。中心段速度剖面与 `u(y)=Δp·y(H−y)/(2 nu L)` 的比较结果：

| 指标 | 结果 |
|---|---:|
| 解析最大速度 | 1.0416667 |
| 模拟最大速度相对误差 | −9.3×10⁻⁹ |
| 速度剖面归一化 L2 误差 | 0.0625% |
| 单位跨度流量误差 | 0.125% |
| 最后写出间隔剖面变化 | 4.8×10⁻⁹ |
| 累计连续性误差 | −5.61×10⁻⁷ |

这验证了单相压力—动量离散能复现通道流，未验证气液界面、润湿或骨架阻碍。

## 静滴表面张力基准

二维圆柱液滴的名义半径为 `R=0.4`，液相/气相均为 `rho=1`，`sigma=1`，`g=0`。解析压差为 `Δp=σ/R_eq`，`R_eq` 根据液相守恒面积计算。数值压差用 `alpha>0.99` 液相和 `alpha<0.01` 气相的体相平均压力差估算。

测试了两种求解器、两种初始化和网格加密：中心点二值初始化；8×8 子单元采样得到的面积分数初始化。面积分数初始化将初始面积误差降至 0.008% 以下。

| 求解器 / 网格 / 初始化 | `nu` / 终止时间 | 液相面积漂移 | 压差误差 | 最终最大速度 |
|---|---:|---:|---:|---:|
| `interFoam` 80×80，二值 | 0.01 / 0.2 | `5.0×10⁻¹³` | −15.30% | 0.2816 |
| `interFoam` 160×160，二值 | 0.01 / 0.2 | `−1.5×10⁻¹¹` | −13.16% | 0.2979 |
| `interIsoFoam` 80×80，二值 | 0.01 / 0.2 | `−8.8×10⁻¹²` | −16.05% | 0.4109 |
| `interFoam` 80×80，8×8 面积分数 | 0.01 / 0.2 | `−4.9×10⁻¹¹` | −12.21% | 0.2721 |
| `interFoam` 160×160，8×8 面积分数 | 0.01 / 0.2 | `−9.9×10⁻¹¹` | −12.58% | 0.3055 |
| `interIsoFoam` 80×80，8×8 面积分数 | 0.01 / 0.2 | `−4.7×10⁻¹¹` | −18.91% | 0.3477 |
| `interFoam` 80×80，`nu=1`，8×8 面积分数 | 1 / 0.5 | `−2.5×10⁻¹¹` | −18.47% | 0.00731 |
| `interFoam` 80×80，`nu=1`，8×8 面积分数 | 1 / 2 | `3.8×10⁻¹²` | −16.07% | 0.00565 |
| `interFoam` 80×80，`nu=1`，8×8 面积分数 | 1 / 10 | `−3.0×10⁻¹²` | −12.56% | 0.00425 |

液相面积和累计连续性误差都很小，但压力跳跃仍比解析值低约 12–22%。二倍网格加密、面积分数初始化以及换用 `interIsoFoam` 都未消除偏差。高黏度 `t=10` 时最大速度降到 0.00425，压力跳跃仍从 `t≈0.02` 的 2.000 上升到 `t≈10` 的 2.186，低于解析值 2.500。仅凭速度减小不能判定静滴平衡通过。

因此，单相通道基准通过，VOF 两相求解可运行且质量守恒良好，但曲率/表面张力平衡测试**未通过**。在改进曲率离散和压力跳跃评估、并完成网格与时间收敛前，不应将其用于 Ag 熔池自由表面、喷溅或熔坑预测。各方案对照见 `static_drop_validation_comparison.png`，数据汇总见 `mesh_sensitivity.csv`，各次完整时间序列保存在 `static_drop_metrics_*.json`。

## 与 Ag 熔池模型的关系

这些测试未包含 Ag 真实密度、黏度、表面张力，未包含 Ti₃SiC₂ 骨架、焓法传热、Ag/骨架润湿、电弧载荷、蒸发或蒸气反冲。64×64×50 熔池数据只有温度、液相分数和骨架掩膜，不能单独确定流动边界条件。即使静滴基准通过，仍须建立包含气相的更大计算域，并对液态 Ag 物性、界面润湿和电弧边界载荷进行文献或实验校准；喷溅预测还需气液耦合，蒸发喷射则需质量损失和蒸气反冲模型。

## 复现

容器镜像为 `opencfd/openfoam-default:2406`，镜像 digest：`sha256:dd5aa20630a55722663bf83ba0cb74870cba130081303e32e3865007fa2aa35a`。工作目录为本报告所在目录。

```sh
python prepare_poiseuille.py
./run_poiseuille.sh
python analyze_poiseuille.py

# 基础网格与二倍网格
python prepare_case.py --resolution 80
./run_case.sh static_drop_case interFoam
python analyze_case.py --resolution 80
python prepare_case.py --resolution 160
./run_case.sh static_drop_case_160 interFoam
python analyze_case.py --resolution 160

# 使用 8×8 子单元面积分数初始化
python prepare_case.py --resolution 80 --case-name static_drop_foam_sub8_80 --solver interFoam --subcell-samples 8
./run_case.sh static_drop_foam_sub8_80 interFoam
python analyze_case.py --resolution 80 --case-name static_drop_foam_sub8_80 --output-tag foamSub8_80

# 可将 --solver 换为 interIsoFoam；高黏度长时间算例可设置 --nu 1、--end-time 10
python prepare_case.py --resolution 80 --case-name static_drop_foam_sub8_nu1_t10_80 --solver interFoam --subcell-samples 8 --nu 1 --end-time 10 --max-dt 0.002
./run_case.sh static_drop_foam_sub8_nu1_t10_80 interFoam
python analyze_case.py --resolution 80 --case-name static_drop_foam_sub8_nu1_t10_80 --output-tag foamSub8_nu1_t10_80

python summarize_vof_tests.py
```

`poiseuille_metrics.json` 和 `poiseuille_profile.png` 保存通道流对照；`static_drop_metrics_*.json` 保存静滴指标及各写出时刻数据；`static_drop_validation_comparison.png` 和 `mesh_sensitivity.csv` 汇总全部十个静滴算例。完整 OpenFOAM 时间步目录可用上述命令重建。
