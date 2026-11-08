# TTMD 熔池候选记录：200 × 200 × 150 Å 细化骨架模型，1.55 ps

本记录对应的不是 917,378 原子的单温 MTP 试算，而是新建的 **TTMD 双温分子动力学**细化骨架模型及其 Round-5 探索轨迹。模型物质区为 200 × 200 × 150 Å，z 向模拟盒延伸至 250 Å，边界为 `p p f`；共有 352,637 个原子，其中 Ag 81,875 个、Ti 135,381 个、Si 45,127 个、C 90,254 个，Ag 质量分数约 50%。固定层为 z=0–20 Å，热浴层 20–40 Å，Newton 分析层 40–150 Å，150–200 Å 与 200–250 Å 是初始为空的溅射/蒸发监测区。

## TTMD 加热方法

这条轨迹使用 Round-5 力拟合开发版 MTP（能量和维里未标定），MD 时间步为 0.1 fs。电子温度用 Ag 占据率加权的 20 × 20 × 15 网格表示，网格间距为 10 Å；XY 周期，z 向电子热流为零。一个中心高斯表面热源先输入 Ag 电子子系统：

\[
q(r_{xy},d)=rac{I_0}{\delta}\exp(-r_{xy}^2/R^2)\exp(-d/\delta),
\]

其中中心为 (100,100) Å，表面 z=150 Å，`I0=5 eV/(ps·Å²)`，`R=30 Å`，电子能量沉积深度 `δ=20 Å`。热源按 Ag 体积映射到 TTM 网格，总积分功率为 3,882.6785 eV/ps；加热到 1.55 ps 时累计输入电子系统 6,018.1517 eV。电子热传导与电子—晶格耦合由 TTM 网格算子处理，再把电子向晶格传递的能量施加给 Newton 区 Ag 原子的动能；Ti₃SiC₂ 没有单独的电子热源，它通过 Ag–骨架相互作用间接升温。因此，这不是单温方案里的“对两相原子直接分配同一份能量”。

这次统计在同一批 1,787 个中心热源区 Newton-layer Ag 原子上做 1 ps 滚动识别，结构每 0.05 ps 采样。联合候选规则为：1 ps 平均 Q6（4.0 Å 邻居截止）<0.4081069；PTM 非 FCC 帧比例≥0.8；去除 Newton 层整体漂移后的 1 ps 解缠绕 MSD\(_i\)≥0.5171679 Å²。候选 Ag 以 3.5 Å 距离连通，XY 周期、z 开放。

## 1.55 ps 结果

同一轨迹的窗口序列从 1.15 ps 到 1.55 ps：1.15 ps 有 424 个液态候选 Ag、89 个连通团簇、最大团簇 97 个原子；1.55 ps 有 686 个候选 Ag、75 个团簇、最大团簇 448 个，占候选 Ag 的 65.3%。最大团簇的几何包络约为 66.2 × 69.4 × 25.1 Å；其中 291/448 个 Ag 距 Ti₃SiC₂ 骨架不超过 3.5 Å。时间序列并非单调长大，最大团簇中间有涨落；因此合适结论是：**TTMD 轨迹确实出现了按同一联合判据识别的、与骨架邻近的连通液态 Ag 候选熔池区**。这比“只有无序原子”更具体，但仍是 Round-5 势下的探索性候选结果，不能作为真实熔池尺寸或抗熔焊性能的定量预测。

1.55 ps 后关闭外部热源并做过关热延续，但到 2.2686 ps Newton 区仍约 5.9 kK、最大原子力达到数百 eV/Å；这段没有形成可信的冷却/凝固证据，也不能据此报告凝固 Ag 桥或熔焊强度。Round-5 势的能量与维里未标定，Ag–Ti₃SiC₂ 热态界面可靠性仍是关键限制。依当前研究顺序，暂不启动新的 TTMD 加热或双触头计算，先继续完成界面 DFT 收敛与分离路径验证。

## 可核查材料

- [1.55 ps 三维联合分类图](runs/ttm_mtp_round5_joint_classifier_to_1p65ps_20261107/joint_melt_1p55ps_3d.png)
- [1.55 ps 联合分类完整 JSON](runs/ttm_mtp_round5_joint_classifier_to_1p65ps_20261107/joint_melt_rolling_report_1p55ps.json)
- [1.15–1.55 ps 团簇演化时间序列](runs/ttm_mtp_round5_joint_classifier_to_1p65ps_20261107/joint_melt_rolling_timeseries_to_1p55ps.csv)
- [加热轨迹的 1.55 ps 状态与 TTM 能量摘要](runs/ttm_mtp_round5_cooldown_from_1p55ps_20261107/summary_restart_1p55ps.json)
- [详细熔池、关热与势函数限制报告](runs/ttm_mtp_round5_cooldown_from_1p55ps_20261107/CONTACT_WELDING_COOLDOWN_REPORT.md)
- [初始 TTM-MD 驱动](runs/ttm_mtp_round5_dt0p1fs_0p05ps_20261107/run_ttm_mtp.py)、[轨迹续算驱动](runs/ttm_mtp_round5_joint_classifier_to_1p65ps_20261107/continue_ttm_mtp.py)、[联合分类分析程序](runs/ttm_mtp_round5_joint_classifier_to_1p65ps_20261107/analyze_joint_melt_rolling_1p55ps.py)

这次记录分支没有打包 352,637 原子结构、Round-5 势文件、逐帧原子/电子轨迹或重启点；它保存 TTMD 方案、分析准则、汇总结果、代码快照与关键图表，避免把摘要误当成完整可复跑归档。
