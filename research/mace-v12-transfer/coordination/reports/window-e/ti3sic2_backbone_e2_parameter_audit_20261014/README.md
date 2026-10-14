# Ti₃SiC₂ 骨架势 E2 归档审计与参数敏感性

**状态：已完成只读 DFT 归档核验与静态参数族敏感性；未拟合候选势。**

仓库没有单独以“Stage E2”命名的根因审计目录。已有 135 原子组件归因与局部环境审计已完成：TSCBOP/DFT 原子映射和实现通过，但 48 原子体相 Ti/C 力失配已存在，114 原子界面派生骨架力 RMSE 为 2.5624 eV/Å。它们没有做参数敏感性。本目录补齐了该项，故根因工作在“骨架项、空间环境、参数族影响”层面完成；原始论文约束、行级参数归属和独立泛化仍未验证。

本目录所有记录均为现有归档、输入 SHA 和输出 SHA 的核对；参数敏感性只用临时副本运行 LAMMPS `run 0`（零步，无时间积分）。原始参数文件哈希保持 `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`。未启动 DFT、MD、TTM、MACE 训练或结构优化，未更改既有 DFT 数据、角色或势文件。

文件：

- `reusable_dft_data_inventory.csv`：可复用 DFT 结果、输入/输出 SHA、单位、方法、收敛和角色。
- `data_quality_and_independence.md`：数据质量、独立性和 stress 缺口。
- `TSCBOP_parameter_sensitivity.csv` / `TSCBOP_parameter_sensitivity_raw.json`：13 个连续参数列的 ±1% 全列方向探针。
- `parameter_sensitivity_analysis.md`：方法、主要结果与限制。
- `constrained_fit_parameter_recommendation.md`：有限的候选参数族建议。
- `current_data_gaps.md` 与 `minimum_additional_DFT_plan.md`：缺口及最小补充计划。
- `role_freeze.json` 与 `fit_start_decision.md`：当前角色冻结和拟合决策。
- `archive_validation.json`：重复执行的哈希验收记录。

本次仅保存本地独立审计目录；没有提交或推送 Git。
