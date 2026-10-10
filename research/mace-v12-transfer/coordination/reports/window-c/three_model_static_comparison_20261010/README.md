# 三模型静态对比交付

本目录保存以 `e455b71f6ef91a0955495001765e659aaee12174` 审计为基线的 C 端静态势函数比较。两种本地模型使用相同公开、已验证参考完成冻结 CPU 推理；MACE-MATPES-PBE-0 许可/元数据受阻，单独记录且没有推理结果。

- `REPORT.md`：主要结果、局限和适用范围。
- `MODEL_SOURCES.md`：官方来源、许可证、权重哈希、元素表及兼容性。
- `V19_REPAIR_PLAN_REVISION_20261010.md`：追加版 V19 决策方案，不覆盖先前的设计。
- `structure_metrics.csv`、`per_element_errors.csv`、`per_atom_errors.csv`：结构、元素与逐原子力误差。
- `ag_response.csv`、`ag_elastic_response.csv`、`backbone_interface_responses.csv`：Ag 体相及 Ti/C 骨架响应。
- `interface_separation_and_layer_forces.csv`：Ag–X 接触分离力与 Ag 层广义 z 力。
- `old_hybrid_vs_dft.csv`：单独抄录的旧混合势历史数据，非本轮计算。
- `reference_inventory.csv`、`model_inventory.csv`、`reproduction_checks.csv`、`provenance.json`：输入、模型和可复现来源记录。
- `compare_static.py`：仅进行冻结推理与输出指标的脚本，不含训练、DFT、MD、LAMMPS 或 TTM 调用。
- `SHA256SUMS.txt`：本交付目录逐文件摘要。

从仓库根目录、使用 MACE 0.3.16 兼容环境重跑：

```bash
python research/mace-v12-transfer/coordination/reports/window-c/three_model_static_comparison_20261010/compare_static.py
```

本次所有结构均按开发诊断或数值控制角色保留；训练资格未变更。封存测试标签、原始模型和 A/B 队列未修改。熔化、凝固及 TTM 适用性未由本静态比较验证。
