# 来源与复现记录

本审计的输入仓库快照为 `linkengde/-` 的提交 `1c4dcc896356860b4d7f678cf9c735e7bbb04bbc`（2026-10-10）。CPU 推理运行时：Python 3.12.14、PyTorch 2.5.1+cpu、mace-torch 0.3.16、ASE 3.29.0、NumPy 1.26.4、SciPy 1.17.1；推理为 float64、CPU。

冻结模型的仓库相对路径及 SHA-256：

- V18：`research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/models/MACE_periodic_v18_clean_core.model` — `757586671174db85acbf7971ccfe92de19cc58378ee2431a797efa6d6163cb62`
- Foundation：`research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model` — `2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012`

逐个 DFT 结构的来源/结果 SHA-256、收敛和 verifier 状态列在 `dft_coverage.csv`；用于六个界面对照的来源与结果哈希列在 `same_geometry_old_potential_comparison.csv`。审计读取固定公开来源，不包含封存验证标签。训练帧清单、云端绝对路径和环境日志不在此发布包中。

复现脚本将所有派生结果写入 `CLOUD_C_AUDIT_OUTPUT` 指定目录（默认系统临时目录），不会覆盖仓库中的已发布 CSV。精确复现应使用上述输入提交中的结构与数据文件。复现可能依赖报告打包后才完成的 A/B DFT 标签；若不存在对应归档，脚本会按报告记录将配对标记为不完整。
