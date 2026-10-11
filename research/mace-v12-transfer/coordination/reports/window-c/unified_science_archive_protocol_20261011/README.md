# A–H 统一科研归档协议与数据索引

本目录落实用户授权的公开科研结果归档规则，并索引 A–H 可从 `linkengde/-` 获取的规范任务资料。现有输入与结果已在规范路径提交时不重复复制；此 C 分支作为读取入口。每个节点的工作产物仍须只提交到自己的分支。

## 文件

- [UNIFIED_ARCHIVE_AND_RECOVERY_PROTOCOL.md](UNIFIED_ARCHIVE_AND_RECOVERY_PROTOCOL.md)：不可覆盖归档、分支规则、SHA-256、GitHub fresh fetch/恢复核验及禁传内容。
- [NODE_DATA_CATALOG.csv](NODE_DATA_CATALOG.csv)：逐节点现有输入/结果路径、完整性清单、来源快照和缺失项。
- [NODE_NOTICES_AND_STATUS.md](NODE_NOTICES_AND_STATUS.md)：A–H 可从仓库读取的正式通知及当前状态快照。
- [CONTROL_EXPERIMENT_MATRIX.md](CONTROL_EXPERIMENT_MATRIX.md)：按现有已批准任务整理的 A/B 界面与骨架对照，以及 D/E/F 3 密度 × 2 seed 液态 Ag 对照。
- [SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json)：分支/来源快照、路径与已有 SHA-256 清单检查结果及未就绪数据。
- [SHA256SUMS.txt](SHA256SUMS.txt)：本目录新增交付文件的校验值。

## 数据传递结论

主分支中的 A/B 已发布结果、C 的门控与执行包、以及 C 分支中的 D–H 任务交接文件均已能从本分支按规范路径读取。数据目录及其原有 `SHA256SUMS.txt` 经核验；本次新增索引不会把已有大目录重复复制或改变来源。

以下项目没有可安全提交的现成数据：Ag.eam.alloy 参数文件（公开再分发状态未核实）、六个 Ag64 液态坐标/轨迹（尚未生成）、G/H 候选 DFT 坐标（尚未生成）、E 独立合成恢复回执（尚未提交）、封存 TEST 数据及未获再分发许可的论文 PDF。缺失项保留为明确阻塞，不用替代文件填充。

本次只新增 C 的协议、索引、节点通知与对照矩阵，不修改 A/B 结果、节点 task card 或 checkpoint；未启动结构生成、MD、DFT、MACE 训练、LAMMPS 或 TTM。D/E/F 暖机仍 HOLD，G/H 仍待命。
