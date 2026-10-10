# 来源与完整性记录

## 基线

- 统一势函数审计提交：`e455b71f6ef91a0955495001765e659aaee12174`。
- 审计目录：`research/mace-v12-transfer/coordination/reports/window-c/v18_unified_potential_audit_20261010/`。
- 审计报告的计算来源快照及推理环境详见其 `PROVENANCE.md`；其中 V18 SHA-256 为 `757586671174db85acbf7971ccfe92de19cc58378ee2431a797efa6d6163cb62`，MACE-MP-0b3-medium SHA-256 为 `2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012`。
- 方案所用的纯相来源限定为 `research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/calculations/` 下的公开 Ag 和 Ti₃SiC₂ 归档；清单每行记录来源结构、结果及摘要哈希、SCF 和 verifier 状态、DFT 口径及数据角色。

## 读取边界

只使用审计报告及 CSV、V18 训练脚本/选择记录/数据 manifest、纯相归档 `summary.json` 与 `verification.json` 元数据，以及公开的 A/B 进度 JSON。未读取 sealed、blind、holdout 或 test 结构/标签；未复制 DFT 能量或力标签。生成 CSV 时只登记来源哈希和数值精度元数据，未将文件系统绝对路径、密钥、训练帧内容或结果标签写入交付物。

纯相归档中的 `energy (force_consistent=False)` 与 `free_energy` 都保留在源结果，但本目录只记录 key 名称和口径说明、不记录能量值。全部 25 条记录的 `dataset_role` 为 `numerical_pure_phase_control_not_training_or_independent_validation`，不得据 SCF 收敛或校验通过推断训练资格。

## 设计依据文件

- `coordination/reports/window-c/v18_unified_potential_audit_20261010/REPORT.md`
- `coordination/reports/window-c/v18_unified_potential_audit_20261010/PROVENANCE.md`
- `coordination/reports/window-c/v18_unified_potential_audit_20261010/per_atom.csv`
- `coordination/reports/window-c/v18_unified_potential_audit_20261010/per_element.csv`
- `coordination/reports/window-c/v18_unified_potential_audit_20261010/ag_response.csv`
- `coordination/reports/window-c/v18_unified_potential_audit_20261010/dft_coverage.csv`
- `periodic_interface_v4/mace_periodic_v18_clean_core/run_v18_training.sh`
- `periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/selection_record.json`
- `periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/training_summary.json`
- `periodic_interface_v4/mace_periodic_v18_clean_core/data/dataset_manifest.json`
- `coordination/progress/window-a.json` and `coordination/progress/window-b.json` as present at the audit baseline.

关键 V18 训练输入的 SHA-256：`run_v18_training.sh` `adf239b551a38191b0017f5bc7e8f3f7de7e18cc40d1903c27002a629b75875c`；`selection_record.json` `d4245492e7d79659a39377cc369f5b43123d876b553eff13043d5e7a534adfe4`；`training_summary.json` `f812f1dff0de3aa4d81e778a860385fd47e78f4936745a075b24913d37b32a8e`；`dataset_manifest.json` `36fb89632293cd5fcd494a4ad7665aeabb575ee1f6cc7633fa913ff7a1f1bdec`。审计结果表哈希可在基线审计目录的 `SHA256SUMS.txt` 核对。基线 `e455b71` 中 A/B 进度 JSON 的 SHA-256 分别为 `6aee3e0d8d9b80acd4aa9da137451b23ca76307828ca27fdfde4de07b4608f30` 和 `536c8c72aceccdbc0bc769439995348ac5bc627a1b35d9ccc08854b88dc6169a`。

本报告的进度状态是 audit baseline `e455b71` 时的快照。后续三模型比较须另核对 `main` 上实际完成且通过验证的归档，不得将 queued/running 项计入有效参考。
