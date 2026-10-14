# Stage E3：Ti₃SiC₂ 骨架势受约束拟合准备

审查基准：`cloud-e-b1-dft` @ `12b42bfbbe0dfc6af936ffe2d1946ed3e27b7a84`。

本目录只包含 E3 审查与方案，不含新的 DFT、LAMMPS 参数扰动、参数优化或模拟结果。E2 已完成的 13 列全列 ±1% `run 0` 敏感性没有重跑。所有现存角色继续冻结为空 TRAIN/DEV/TEST。

## 结论

- E2 的 13 个归档文件均存在，清单中的 12 个载荷 SHA-256 全部通过。
- TSCBOP 文件 SHA-256 为 `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`；未修改。
- 文件可见的参数共享模式已盘点，但 Plummer–Tucker 原文/附录不在本地；外部 27 条记录、378 个数值字段的比对视为外部已核验，E 本地未独立核验。
- E2 扰动记录只有结构/元素汇总 RMSE，没有逐原子参数响应向量。不能据此构造逐行 Jacobian、矩阵秩或奇异值，也不能确认可辨识参数。
- 当前不建议任何实际试拟合：`START_TRIAL_FIT=NO`。先经批准补做 tied-block 静态 `run 0` 力响应，再补齐关键 DFT 数据并冻结数据角色。

## 文件

- `E2_REVALIDATION_AND_SOURCE_STATUS.md`：E2 哈希、基线和来源核查。
- `TERSOFF_PARAMETER_MAPPING_AND_BLOCKS.md`、`TERSOFF_PARAMETER_BLOCK_MAP.csv`：参数函数角色及文件中可观察到的共享模式。
- `PARAMETER_RESPONSE_LIMITS_AND_RUN0_PLAN.md`：可辨识性边界及待批准的最小响应测试。
- `FIT_OBJECTIVE_AND_DATA_ROLES.md`：目标函数和无泄漏角色方案。
- `MINIMUM_DFT_PLAN.csv`：去重后的补充静态单点建议。
- `FIT_DECISION.md`：决策及对用户问题的逐项回答。
- `E2_FILE_HASHES.csv`：13 个 E2 包文件的本地 SHA-256 记录。
- `SHA256SUMS.txt`：本目录交付文件清单。
