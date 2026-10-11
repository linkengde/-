# A–H 节点通知与当前状态

日期：2026-10-11
通知来源：云端 C
通知落点：C 分支 `cloud-c-three-paper-crystallography-interface-audit-20261011` 的统一归档协议目录

> 本文件为供各节点从仓库读取的正式通知。它不修改节点分支、task card、checkpoint 或任务队列。节点可在自己的 checkout 只读执行 `git fetch origin cloud-c-three-paper-crystallography-interface-audit-20261011`，再用 `git show FETCH_HEAD:<catalog 中的路径>` 读取所需文件；如需落盘，应写入全新的节点专属 staging 目录，不合并 C 分支、不覆盖已有文件。未收到节点回执前，不表示其已阅读或确认。

## 全节点共同执行要求

请先阅读同目录 `UNIFIED_ARCHIVE_AND_RECOVERY_PROTOCOL.md`、`NODE_DATA_CATALOG.csv` 和 `CONTROL_EXPERIMENT_MATRIX.md`，再按 catalog 中的规范仓库路径取数。只处理自己节点和自己授权 task。首次归档和每阶段结果必须有唯一目录名、README、MANIFEST、完整逐文件 SHA-256、正常提交/推送记录，以及 GitHub fresh fetch 到新目录后的恢复校验。结果文件、归档状态及数据许可必须如实标识。

此 C 分支包含已发布的 A/B 结果、C 任务包及 G/H 交接材料；仓库原有规范路径即为共享数据源，不要为“共享”再复制整棵目录或重写历史。不同节点的工作仍只发布到本节点自己的分支。需要的内容若未列为已存在，先向 C 报缺项，不得自行生成或猜造输入。

## 当前节点状态表

| 节点 | 仓库中的数据/最近证据 | 当前状态与门槛 |
|---|---|---|
| A | Ag–Si parent 与 ±0.20 Å 三点在 main 快照 e91ef26 中 | 既有任务照常；本轮不改任务或 checkpoint |
| B | Ti/C 四个位移最终结果在 34b8be9 归档 | 已完成；数值控制角色不变 |
| C | C 协调分支含统一归档协议及 D 首项暖机签署 | 汇总节点回执；只释放 D 首项，不代跑节点任务 |
| D | owner 分支 `2c8b6f41ec8a7605ddc06f7218014c55957166ac`；本机势 SHA、ASE 解析、资源及归档回执已核验 | C 已签署唯一首项的 1000 步暖机；D 必须以签署 SHA 重跑实时 preflight 并 PASS 后才能启动 |
| E | owner 分支 `39734568216e5a3d9853d7d052d28f9799172598`；势 SHA、ASE 元数据解析及独立合成远端往返回执已核验 | 环境证据项通过；node preflight 和科学任务仍 HOLD。cutoff 跨半盒长由已审查的 ASE 多周期像邻居表支持 |
| F | owner 分支 `2ddb39a899decd6b3cc8d6c002d48291c9a5e632`；势 SHA、ASE 元数据解析、资源身份预检及归档回执已核验 | 环境子检查通过；F host release=false，科学任务仍 HOLD |
| G | 本轮远端分支头未核实；HANDOFF 已在 C 分支 | 待命；无结构、无新增任务 |
| H | 本轮远端分支头未核实；HANDOFF 已在 C 分支 | 待命；无结构、无新增任务 |

## 各节点通知

- **A：** 你已发布的结果和现有任务状态原样保留。继续使用 main 中现有规范结构/参考路径；后续经批准的结果使用 A 自己的分支和不可覆盖归档。不要改写旧归档，也不要重新启动已完成工作。
- **B：** Ti/C 四个对称位移结果及报告（来源提交 `34b8be9`）保留原样，按 catalog 的路径读取。它们的数据角色不因新归档协议改变；不得自行加入训练。
- **C：** 汇总 A–H 的归档回执和数据清单；每次核对节点分支当前头、归档外层/内部 SHA、恢复结果和科学状态。C 分支提交只用于协调文件，不代替节点输出归档。
- **D：** C 已签署仅限 `C-D-AGL-1250K-rho097-seed01` 的最多 1000 步、0.5 fs 暖机，配置 SHA 为 `4aee9d64859736f35459259625598f767d5f3e0d8b33bc71f0abb5baf31639e2`。D 必须先在当前会话按 `codex_ephemeral_ag64_warmup_release_20261011/C_RELEASE_ATTESTATION_D_WARMUP_20261011.md` 重跑只读预检并取得 PASS；执行器会再次检查并在完成后归档/恢复核验。不得自动续算或启动 D 第二项。
- **E：** 使用 C 分支中的 `START_E.md` 与门控脚本；仅在自己的 E 分支提交。E 的独立 synthetic fixture 已完成 fresh fetch、恢复和 SHA 核验（owner 提交 `39734568216e5a3d9853d7d052d28f9799172598`）。本机势 SHA 和 ASE 元数据解析也通过；旧配置下的预检为 HOLD。cutoff 相对半盒长的检查已有 C 的 ASE `PrimitiveNeighborList` 多周期像审查依据；任务本身仍 HOLD，因当前 C 配置不释放 E。
- **F：** 使用 C 分支中的 `START_F.md` 与门控脚本；仅在自己的 F 分支提交。F 的两项当前预检及报告远端读回已经核验（owner 提交 `2ddb39a899decd6b3cc8d6c002d48291c9a5e632`），势 SHA、解析、资源/路径及 archive 子检查通过；F 的 `C_host_release=false` 且 task 不在放行范围内。不得启动 F 两项任务。
- **G：** 仅取得 `HANDOFF_G.md` 和任务矩阵供后续准备；当前没有结构输入和新任务放行。保持待命，不生成结构、不启动 DFT。
- **H：** 仅取得 `HANDOFF_H.md` 和任务矩阵供后续准备；当前没有结构输入和新任务放行。保持待命，不生成结构、不启动 DFT。

## 状态快照边界

本次只读查询的 GitHub 头：`main=e91ef26`；A 报告引用该 main 快照，B 最终结果来源 `34b8be9`。D/E/F 的 owner 分支最近核验头分别为 `2c8b6f41ec8a7605ddc06f7218014c55957166ac`、`39734568216e5a3d9853d7d052d28f9799172598`、`2ddb39a899decd6b3cc8d6c002d48291c9a5e632`。D/E/F 详细验收结论见 `NODE_EVIDENCE_UPDATE_E_F_20261011.md` 与 `D_WARMUP_RELEASE_SUMMARY_20261011.md`。节点启动任何获批工作前仍须在自己的会话刷新分支头和资源。

本次查到的远端分支列表未给出 G/H owner 分支头，故记作 `UNVERIFIED`。没有通过本通知建立或推送 G/H 分支。

归档权限不等于科学执行授权。D/E/F Ag64 暖机逐项受独立 gate。cutoff 大于半盒长已由冻结的 ASE 多周期像邻居表实现支持；该静态代码审查不代表 64 原子下的熔点、扩散或凝固物性验证。G/H 待命，A/B 的既有状态不变。不同对照轴及 D/E/F 的既有 3 密度 × 2 seed 分工见 `CONTROL_EXPERIMENT_MATRIX.md`；该矩阵不新增 task，也不解除放行门。
