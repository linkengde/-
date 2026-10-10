# 分支关系和安全登记提案

handoff bb8ad7e4b725b634de2e2198d141df2bfe9d78b9 与 main 684e9a887547404d959e2a7b679c9596ea134a8a 共同祖先 bcb51e5e05493e846ef3bf2cc86081dba8b11abe；handoff有1个C提交，main有3个A/B进度提交。交接尚未进入main。

本轮C专属分支 cloud-c-preflight-pilot-readiness-20261011 从原C交接分支建立，用普通merge 5106470ab47f3ea8b8fe9dd5c82d0bd04c0fe471 合入审查main，无冲突。A/B公开文件内容与该main一致；运行工作区未动。新科学材料只增加本报告目录，不改旧交接。只正常push新C branch，不push/force-push main。

D提供的f9544cd00046e7fd4be75fec4d6fe13e50e8cdab是仓库快照/历史定位信息，不是D资源报告SHA；该commit实际是B进度提交，不能据此声称D报告已归档。G/H本地bcb快照较C审查main旧，可只读fetch再核对自身branch；运行input/checkpoint保持原位，不checkout覆写活动输入。

后续安全登记/合并提案（本轮不执行）：
1. review旧handoff、再review本次C新报告；本C分支差异含旧17份交接与本目录，两个C报告提交分别可查。
2. 每个owner采用自己的cloud-x-b1-dft分支、task卡/manifest文件；只登记自己6个task ID。C在自己汇总文件维护主表，不写A/B任务卡。
3. 状态按PLANNED→SOURCE_VERIFIED→GENERATION_AUTHORIZED→STRUCTURE_QA_PASS→DFT_PILOT_AUTHORIZED→NUMERICAL_PASS→FINAL_LABELS_AUTHORIZED→ARCHIVE_VERIFIED，独立user approval不跳过。
4. 六个pilot是原30任务的子集；数值mesh control使用task_id/mesh_id作为run身份，不加训练父族/正式构型数量。
5. 提交/合并前fetch并检查目标paths/branch；同文件并发更新、非快进或冲突即停止该提交并报告。不强推、不删除其他cloud结果、不用ours/theirs丢弃数据。只增加各自文件后正常review合并main，不让多个host直接改同一汇总文件。
6. main上A/B进度增长只作为新公开快照；本报告来源固定以上SHA，不声称实时状态。保持所有活动任务不变。

旧交接术语纠正：TSC 10×10×4基准是已审计48atom reference cell，不是primitive cell。在新报告更正，不覆盖原HANDOFF。

## 发布前安全同步记录

发布前fetch发现main已推进至284285cec6a7a5691ffc01fd9616d7760dc4dac4，新增B的公开归档/状态提交d359f702和284285c，未改C报告路径。C用普通merge d868badbee761afd3566c2d1c3fac150968f2e65 安全纳入，零冲突，保持这些公开文件与main一致；没有重写main或修改B工作区/数据角色。本轮不展开新增标签分析，原报告科学/资源证据快照仍为684e9a8加owner反馈。

最新main相对共同祖先bcb51e5领先5个A/B提交；旧handoff bb8仍不属于main历史。本C分支包含最新已审查main和原handoff，main上的其他人结果完整保留。
