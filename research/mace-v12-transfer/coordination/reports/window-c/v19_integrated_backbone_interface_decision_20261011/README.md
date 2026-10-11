# 云端 C：骨架—Ag–Si界面综合修复决策

本目录是2026-10-11独立追加报告，不覆盖任何旧审计、原数据、模型或A/B记录。发布在C专属分支`cloud-c-v19-integrated-decision-20261011`。

1. 先读[REPORT.md](REPORT.md)：四类受力证据和结论。
2. 读[V19_REPAIR_DECISION.md](V19_REPAIR_DECISION.md)：对照臂、能量口径及候选停止门。
3. 读[DATA_ROLE_AND_INDEPENDENCE.md](DATA_ROLE_AND_INDEPENDENCE.md)及[candidate_data_requirements.csv](candidate_data_requirements.csv)：控制角色不变、未来来源隔离与成本。
4. 读[EXECUTION_ORDER.md](EXECUTION_ORDER.md)：建议顺序，不是执行授权。

复核数据：`interface_per_atom.csv`312行、`interface_per_element.csv`24行、`interface_metrics.csv`6行、`interface_response_comparison.csv`14行；`backbone_*_reused.csv`原字节复用C已发布结果。`pure_phase_role_audit.csv`逐条26个公开控制记录，当前训练资格全部NO。`PROVENANCE.json`及来源review/verification保存固定提交、SHA和工作边界。

`compare_interface_frozen.py`是本轮来源验证与两项缺失冻结推理的复现脚本；本轮已完成，**不要为验收再次运行它**。交付复核使用`validate_delivery.py`（仅数组/哈希/CSV检查，不加载模型）。源数据仅按明确白名单读取，不扫描sealed；不给训练组包入口或DFT/MD运行入口。校验：在本目录执行`sha256sum -c SHA256SUMS.txt`。

`V19_TRAINING_RELEASE=HOLD`、`NEW_SOLID_DFT_RELEASE=HOLD`。六液态pilot继续原有存储/暖机放行门；G/H待命，其余24HOLD。当前文件不是任何科学计算的新授权。
