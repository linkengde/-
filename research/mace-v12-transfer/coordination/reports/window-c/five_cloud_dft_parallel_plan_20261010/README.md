# 五云端 DFT 交接包

此目录为 C 生成的 Ag–Ti₃SiC₂ B1 DFT 计划与五份独立交接文件。所有构型都只是提案；没有坐标或可运行输入，当前不允许启动任何 DFT/MD/训练。

## 文件

- [`HANDOFF_D.md`](HANDOFF_D.md)、[`HANDOFF_E.md`](HANDOFF_E.md)、[`HANDOFF_F.md`](HANDOFF_F.md)、[`HANDOFF_G.md`](HANDOFF_G.md)、[`HANDOFF_H.md`](HANDOFF_H.md)：分别单独复制给相应云端。每份含本机六个 task ID、来源/parent/seed/结构状态、环境只读预检、DFT/k 点要求、分阶段授权、归档与停止条件。
- [`MASTER_TASK_MATRIX.csv`](MASTER_TASK_MATRIX.csv)：全局唯一 30-task 主表，包含每个候选的预计费用、hash 状态、数据角色、k 点和验收门。
- [`RESOURCE_AND_TIME_BUDGET.md`](RESOURCE_AND_TIME_BUDGET.md)：实际资源可见性、建议资源门、pilot 与全批成本估算。
- [`PILOT_RELEASE_GATES.md`](PILOT_RELEASE_GATES.md)：1250 K 液态 Ag 六个 pilot 的 subset 关系、来源/液态/k-point/资源/数据验收门。
- [`DATA_IDENTITY_AND_SPLIT.md`](DATA_IDENTITY_AND_SPLIT.md)：parent-family、seed、SHA 去重和训练/开发/测试角色规则。
- [`five_cloud_task_matrix.csv`](five_cloud_task_matrix.csv) 与 [`resource_preflight.csv`](resource_preflight.csv)：简洁分配表及当前 D–H 可见资源/任务状态。
- [`REPORT.md`](REPORT.md)：总计划、限制、成本、执行阶段。
- [`SHA256SUMS.txt`](SHA256SUMS.txt)：除自身外本目录交付文件的 SHA-256。

## 使用顺序

1. 先由 C/协调者核对 `MASTER_TASK_MATRIX.csv` 和各 HANDOFF，确认 task ID/parent/seed 不重复。
2. 仅将每个云端自己的 HANDOFF 发给对应 owner。D/E/F 报告 pilot 指定 ID；G/H 在 pilot 阶段预留。
3. 云端 owner 先做只读环境和当前任务状态预检。活动任务/checkpoint 保持原样；所有结构/hash 当前仍 `STRUCTURE_NOT_READY`，需报告来源和实际缺口。
4. 补齐 1250 K 液体密度依据、采样器许可及真实液体候选结构后，C 复核具体输入和收敛预算。结构生成 MD/AIMD 与六个 DFT pilot 都须用户分别明确批准。
5. 用户批准 pilot 后，D/E/F 各仅处理自己分配的两个 pilot task ID。C 审计收敛、k 点、液态性、重复、哈希、耗时和资源；未获进一步批准，其他 24 行不得运行。
6. 六点全部检查后，C 报告 pilot 结果和修订预算；只有用户再次批准，才能释放余下 24 项。

`HANDOFF_READY = YES` 只表示交接材料和清单可供只读预检查使用。`DFT_EXECUTION_READY = NO`：没有候选坐标、输入文件或 SHA，密度/采样器仍有阻塞，D–H 实际环境也未由 C 核验。本目录不是运行授权。
