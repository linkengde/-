# Ag–Ti₃SiC₂ DFT修复：聊天界面生成文件直接提交仓库

## 文件究竟从哪里来

- **本次提交的生成位置**：本次ChatGPT聊天界面的工作目录 `/mnt/data/Ag_Ti3SiC2_DFT_Reference_Repair_20261016/`，不是E/C窗口的 `/workspace/`。
- **提交工具/目的地**：ChatGPT当前窗口通过已连接的GitHub仓库接口，提交到仓库 `linkengde/-` 的独立分支 `audit-nep-ace-dft-reference-repair-20261016`，目录即本README所在目录。
- **原始科学数据来源**：仓库旧提交 `499ddf0870235594444abad780582ea4fcc89d61`，包含V16训练集的第0帧和历史同几何DFT资料。**原始DFT不来自ChatGPT生成文件**；此窗口只是读取原始数据并离线统计、重建无标签输入和整理计划。
- **本次写入类型**：`CHATGPT_GENERATED_DERIVED`，数据表和无标签结构均为本聊天生成的派生材料，不能被标为新的DFT结果。
- **执行状态**：NO_NEW_DFT / NO_TRAIN / NO_MD / NO_TTM。DFT标签冲突未解决，NEP/ACE训练数据未批准。

## 清晰表格与含义

| 仓库文件 | 聊天界面源文件 | 类型 | 状态 |
|---|---|---|---|
| KPOINT_LADDER.csv | /mnt/data/Ag_Ti3SiC2_DFT_Reference_Repair_20261016/KPOINT_LADDER.csv | k点网格参数的离线计算表 | 新网格仅已准备，未计算SCF |
| DFT_ACQUISITION_QUEUE.csv | /mnt/data/Ag_Ti3SiC2_DFT_Reference_Repair_20261016/DFT_ACQUISITION_QUEUE.csv | 后续DFT任务规划 | 待审查，不是已运行队列 |
| DFT_SOURCE_REUSE_RECONCILIATION.csv | /mnt/data/Ag_Ti3SiC2_DFT_Reference_Repair_20261016/DFT_SOURCE_REUSE_RECONCILIATION.csv | 原始数据复用/去重统计 | 只覆盖先前已读取的资料范围 |
| Ti_gap_2p3_frozen_32atoms_UNLABELED.extxyz | /mnt/data/Ag_Ti3SiC2_DFT_Reference_Repair_20261016/Ti_gap_2p3_frozen_32atoms_UNLABELED.extxyz | 无能量/力标签冻结结构 | 可供后续重新验证，不是DFT输出 |
| CONFLICT_AND_EXECUTION_GATE.md | 由本聊天生成，依据前轮执行报告/仓库的原始DFT摘要 | 参考冲突审计和执行门槛 | 对照尚未收敛 |

**注意**：以上是“聊天界面生成 -> GitHub独立分支”的文件流。以前E窗口写入仓库的报告不能被表述成聊天界面上传；直接从GitHub读取的DFT结果也不能称作在本窗口新算出的DFT。

## 原始GitHub依据

- `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v16_interface_energy/data/train.extxyz` 第0帧：`Ti_gap_2p3`，Ag10Ti12Si2C8、32原子。
- `research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v17_rank_recovery_parallel/calculations/Ti_gap_2p3_rank_recovery_v17_01/`：同几何Γ点重算，包含summary及verification。
- `research/mace-v12-transfer/coordination/reports/window-c/v19_controlled_repair_design_20261010/pure_phase_source_eligibility.csv`：纯相DFT数值控制来源、角色与SHA。
- `research/mace-v12-transfer/coordination/reports/window-a/v19_Ti3SiC2_bulk_mesh_comparison/`：48原子骨架k网格参考。

旧数据均保持只读；本分支不覆盖或修改它们。
