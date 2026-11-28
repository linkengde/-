# 材料物性来源审计

## 本阶段使用的统一参数

所有纯相解析基准和RVE结果均使用同一组300 K常数热导率：Ag `k=429.00 W/(m·K)`；Ti₃SiC₂ `k=37.25 W/(m·K)`。Ag候选CSV备注为“User-supplied baseline values attributed to NIST/CDC thermophysical table”。Ti₃SiC₂由仓库保存的相关式 `k_W_m_K(T) = 38.6 - 0.0045*T` 在300 K计算，报告拟合区间为25–1200 °C、R²=0.86。本轮等效热导率线性计算只用 `k`，不使用Cp、密度或相变参数。

## 可追溯信息和限制

- 求解器定义对照自 `main` 提交 `fbc7fbd1081fcb4a2a86a671188c3327a577631f`（Blob `708df50bd80f4fe84fbb2b38e16026f97f3da9fd`）以及 `codex/single-temperature-study` 提交 `703b63ccf840dc0f8d0e85349221c8457163417f`（Blob `43826c0aa746659a36c3b4b03de5ad9ff4e6b8b1`）。原有求解器是显式瞬态焓法；本阶段独立使用稳态FVM模块，不修改原求解器或MD代码。
- 仓库把Ag热导率归于NIST/CDC资料，但缺表名、版本、行号/页码、URL和DOI；Ti₃SiC₂拟合归于Barsoum et al. (1999) Sec. 3.4/Fig. 6并由Qiu et al. (2020) Eq. 24复述，但提取文档没有完整书目与DOI。
- Ag–Ti₃SiC₂界面热阻没有可核验实测值。`10⁻¹⁰`、`10⁻⁹`、`10⁻⁸ m²·K/W`只用于算法和灵敏度演示，不是实验物性。

因此，解析基准证明了求解器实现与声明参数一致；绝对有效热导率仍为候选物性条件下的结果。完成完整引文和原始数据核对前，不能将其称为已审计的材料定量预测。
