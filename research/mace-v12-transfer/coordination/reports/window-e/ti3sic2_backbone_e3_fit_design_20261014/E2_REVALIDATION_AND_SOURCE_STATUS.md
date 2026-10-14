# E2 归档复核与论文来源状态

## 审查范围与基准

- 仓库：`linkengde/-`
- 分支：`cloud-e-b1-dft`
- 用户指定基准：`12b42bfbbe0dfc6af936ffe2d1946ed3e27b7a84`（当前 HEAD）
- E2 目录：`coordination/reports/window-e/ti3sic2_backbone_e2_parameter_audit_20261014/`
- E2 归档包含 13 个文件：12 个受 `SHA256SUMS.txt` 校验的载荷及该清单本身。12/12 项通过；无缺失或不匹配。

E2 包含数据清单、质量/角色审查、13 列敏感性 CSV/JSON、参数建议、数据缺口、补充 DFT 计划、角色冻结、审计记录和 README。逐文件 SHA 见本目录 `E2_FILE_HASHES.csv`。没有重新执行 E2 的 13 列 ±1% 敏感性计算。

## 归档记录复核

`archive_validation.json` 记录的原归档检查为：纯相控制包 29 项通过、组件力归因 53 项通过、界面弛豫包 17 项通过；7 个体相基线 mesh 子包和 4 个 Ti/C 位移子包也通过各自 SHA 检查。TSCBOP 原文件 SHA 与预期一致。

E2 的敏感性是 13 个连续数值列分别对全列非零记录统一乘以 0.99/1.01；共 26 个参数副本、156 个结构/副本 `run 0` 评估，6 个结构案例。它是已有结果，本轮只读复核其清单与限制，没有重做计算。

归档基线力 RMSE：48 原子单一父族的 Ti/Si/C 分量约为 2.2206/0.0062/2.1290 eV Å⁻¹，总体约 1.9941 eV Å⁻¹；114 原子界面派生 TSC-only slab 的 Ti/Si/C 约为 2.6242/2.6072/2.4433 eV Å⁻¹，总体约 2.5624 eV Å⁻¹。这里的 48 原子输入是 CIF 派生的固定近父几何，不应称作已证实的 PBE 平衡晶体。

## 文献状态

用户提供的独立审查结论称：Plummer 与 Tucker，*Physical Review B* 100, 214114 (2019)，DOI `10.1103/PhysRevB.100.214114`，附录 Table V–VII 与仓库文件的 27 条记录、378 个数值字段一致。E 在当前工作区的仓库与附件中未找到该论文或附录原文；本轮未能独立读取原文，故状态为：

`SOURCE_STATUS=EXTERNALLY_VERIFIED_PER_SUPPLIED_AUDIT; E_LOCAL_VERIFICATION=NOT_VERIFIED`

原 TSCBOP 文件 SHA-256 为 `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`。E 检查了文件头、27 行及 14 个数值字段的格式；这不是对论文参数的独立逐项比对。

论文/原文件中的参数值只是原始参考点 `p0`。它们不构成拟合允许上下限。E2 的 ±1% 是局部敏感性扰动步长，也不是参数边界。当前没有来源充分的数值拟合上下限。

## 本轮基线状态

`BASELINE_REPRODUCED=PASS`（限于 E2 已归档并复核的基线输出/坐标对应与 RMSE；本轮没有重新运行 LAMMPS）。
`E2_PACKAGE_SHA=PASS`。
`SOURCE_TABLES_E_LOCAL=NOT_VERIFIED`。
