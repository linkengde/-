# 云端 C：Ag–Ti₃SiC₂ 统一势函数审计

本目录归档 V18 骨架力退化、界面误差、旧候选混合势对照、固态 Ag 响应和 DFT 覆盖审计。审计报告给出数值定义、证据限制、液态/凝固数据缺口及后续修复和验证建议。

- 主报告：[`REPORT.md`](REPORT.md)
- 逐原子误差：[`per_atom.csv`](per_atom.csv)
- 分元素误差：[`per_element.csv`](per_element.csv)
- 基础 MACE/V18 固态 Ag 响应：[`ag_response.csv`](ag_response.csv)
- 同几何旧混合势与 MACE 对照：[`same_geometry_old_potential_comparison.csv`](same_geometry_old_potential_comparison.csv)
- DFT 资产覆盖：[`dft_coverage.csv`](dft_coverage.csv)
- 来源、输入模型哈希与复现说明：[`PROVENANCE.md`](PROVENANCE.md)

四个 Python 脚本提供冻结模型推理与路径筛选审计。它们默认把生成文件写入系统临时目录；设置 `CLOUD_C_AUDIT_OUTPUT` 可更改输出位置，设置 `CLOUD_C_REPO_ROOT` 可指定输入仓库根目录。它们不会启动训练、DFT、MD 或 TTM。`compare_models_hybrid.py` 只读取既有候选02静态指标，不运行 LAMMPS。

运行环境要求和精确版本见报告及 `requirements-mace.txt`。本次发布包已排除云端运行日志、本地安装脚本、完整原始审计 JSON，以及含历史 `holdout` 配置字样的训练帧清单；未发布封存标签。

校验归档文件：

```bash
sha256sum -c SHA256SUMS.txt
```
