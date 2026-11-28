# Stage 4 material-property audit status

This is the same candidate-property audit used by `phase_change_validation/`; it is copied here so solver-performance results carry their material-data limits. Stage 4 made no thermophysical parameter changes. A local repository search still found no primary Ag tables or Barsoum source PDF/full citation, and the network literature lookup was blocked by the execution proxy (Crossref request returned HTTP 403). No DOI, source table, or numeric data has been guessed. The NO-GO decisions below therefore remain active.

---

# Ag/Ti₃SiC₂ 热物性来源审计

## 审计结论

本阶段的**数值方法验证通过**，但用于实际工况预测的**材料物性审计未通过**。Ag 数据只有仓库内的候选值及“NIST/CDC 表”归属说明，缺少可定位的表名、版本、行号/页码、原始 URL 和 DOI。Ti₃SiC₂ 关系可追溯到文献中报告的相关式，但仓库没有完整书目和 DOI，也没有原始数字化数据及不确定度。本次执行环境内未发现相应原始 PDF；因此没有补写或猜测 DOI。

本报告把可复核的仓库参数、拟合关系和数值假设分开记录。所有相变结果仅为**候选物性条件下的探索性数值验证**，不是已经完成文献审计的材料预测。

## 参数清单与单位核对

| 参数 | 仓库候选值/模型 | 换算及适用范围 | 数据属性与限制 |
|---|---|---|---|
| Ag 摩尔质量 | 0.1078682 kg/mol | 用于 J/mol→J/kg | 候选输入；此次未独立核对原始表 |
| Ag 熔点 | 1235.08 K | 固液转换温度 | 候选输入；不是由 MD 势函数预测 |
| Ag 熔化潜热 | 11297 J/mol = 104729.6608 J/kg | 除以 0.1078682 kg/mol | 候选输入；单位换算复算一致 |
| Ag 固态密度 | 10490 kg/m³ | 焓体积项使用此值 | 候选值；两相均固定用固态参考密度 |
| Ag 液态密度 | 9346 kg/m³，熔点处 | 仅留作候选记录 | **没有在本焓法中使用**；不模拟膨胀、体积变化或质量迁移 |
| Ag 固态 Cp | 25.4 J/(mol·K) @300 K；31.8 J/(mol·K) @1235.08 K | 235.4725、294.8042 J/(kg·K)；两点间线性插值 | 候选锚点，源文件称归于 NIST/CDC；具体表格和温区不完整 |
| Ag 液态 Cp | 33.47 J/(mol·K) @1235.08 K | 310.2861 J/(kg·K) | 单一锚点；实现中保持常数，不能据此声称温度趋势已测得 |
| Ag 固态热导率 | 429 W/(m·K) @300 K；175 W/(m·K) @熔点 | 两点间线性插值 | 候选锚点；无完整原始出处和适用范围 |
| Ag 液态热导率 | 175 W/(m·K) @熔点 | 实现中保持常数 | 仓库标注为 tentative；没有随温度变化的独立数据 |
| Ti₃SiC₂ 密度 | 4520 kg/m³ | 固定参考密度 | 仓库模型输入；不耦合热膨胀 |
| Ti₃SiC₂ Cp | `Cp = 164.4 − 16419/T` J/(mol·K)，T 用 K | 除以 0.195708 kg/mol 得 J/(kg·K)；解析积分用于焓 | 仓库归于 Barsoum et al. (1999) Sec. 3.3/Fig. 5，Qiu et al. (2020) Eq. 25 复述。原始 DSC 报告 100–1000 °C，原文外推至约1200 °C；实现从300 K至1473.15 K，低端/高端均有超出直接 DSC 区间之处 |
| Ti₃SiC₂ 热导率 | `k = 38.6 − 0.0045 T` W/(m·K) | @300 K 为37.25；@1473.15 K 为31.9708 W/(m·K) | 仓库归于 Barsoum et al. (1999) Sec. 3.4/Fig. 6；报告拟合范围25–1200 °C、R²=0.86。超出该范围不作外推 |
| Ti₃SiC₂ 各向异性 | 本求解器用标量 k(T) | 无张量分量输入 | 仓库没有本阶段可核验的各向异性数据；本模型不表示其各向异性已被排除 |
| Ag/Ti₃SiC₂ 界面热阻 | 主三维相变算例 R=0 | 理想热接触基线 | 实测值缺失；R=0 不代表真实界面必然无热阻 |
| 界面热阻敏感性值 | R=10⁻¹⁰、10⁻⁹、10⁻⁸ m²·K/W | stage2 RVE 的平均对角热导率在10 μm模型相对R=0分别变化−0.217%、−1.539%、−5.742%；200 μm模型为−0.011%、−0.112%、−0.905% | **数值敏感性参数，非实验物性**；结果来自已有稳态验证报告，本次三维相变主算例未加入有限R |

## 来源追踪状态

- Ag 候选 CSV 备注为“User-supplied baseline values attributed to NIST/CDC thermophysical table”，但具体出版物、版本、行页、URL/DOI 均未保存，故不能逐值复核。
- 仓库 JSON 记录 Ti₃SiC₂ 相关式来源为 Barsoum et al. (1999) 第3.3节图5（Cp）和第3.4节图6（k），并注明 Qiu et al. (2020) 方程25/24复述。当前 checkout 中缺完整参考文献条目和 DOI；不把二手综述当作独立实验数据。
- 仓库还提到 Togo et al. (2010)、Wang and Zhou (2004)、Son et al. (2016) 的计算结果作为交叉检查；它们不是高温 Cp/k 的独立实验锚点，也不能证明2000 K可用。
- 上游物性 JSON 中的2000 K外推不被本求解器接受。本阶段 Ti₃SiC₂ 温度守卫设为1473.15 K（1200 °C）；超过即报错。此上限是采用的参数关系范围守卫，不代表材料稳定上限。
- 没有使用 MD/MACE 预测的熔点作为 FVM 熔点。

## GO/NO-GO

| 范围 | 判断 |
|---|---|
| 潜热单位换算及程序中的 kg/J/m³转换 | GO；显式换算复算一致 |
| 候选参数下的单温焓法代码验证 | GO；见 `STEFAN_VALIDATION.md` 与小算例结果 |
| 实际 Ag/Ti₃SiC₂ 定量热响应预测 | **NO-GO**；Ag 表格出处不完整，液态锚点稀疏，界面热阻无实测值 |
| Ti₃SiC₂ 高温各向异性及组织稳定性 | **NO-GO**；当前使用标量热导率，且未建立相分解模型 |

因此，任何熔化深度/凝固时间仅可作为这组输入参数的条件性敏感性结果，不应标成已验证的触头实际工况预测。
