# Ti₃SiC₂ 独立 DEV 父族来源审计

日期：2026-10-20  
分支：`cloud-e-b1-dft`  
基线提交：`bfb09bcc69ab83b4e5147638ec91f2bc981f0a6c`

## 结论

`DEV_PARENT_SOURCE=NOT_FOUND`  
`DEV_DFT_READY=NO`  
`TRAIN_DEV_TEST_ROLES=UNCHANGED`  
`TEST_BOUNDARY=RESPECTED`  
`DFT_EXECUTION_READY=NO`

本轮仅做来源与运行环境审计，没有构建或改动结构，没有启动 DFT、MD、优化或拟合。

## 来源审计

E 批次任务清单 `tsc_backbone_dft_batch_20261017/PUBLIC_TASK_MANIFEST.json` 中，任务 5（DEV 参考点）和任务 6（其确定性位移响应点）仍是 `HOLD_NO_VERIFIED_INDEPENDENT_NONTRAIN_SOURCE`，结构路径和 SHA 均为空。该清单的 SHA-256 为 `9997106fd1b2f8963453f9163ebbb98e4ab420e581d7c4488a4ca39d2c50b6ef`。

仓库 Ti₃SiC₂ 晶体源扫描只找到：

- `research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif`
- SHA-256：`6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3`
- 来源：COD 9009647；记录引用 Kisi 等人 1998 年的晶体结构。

四个已完成的 100 原子表面点均从这份 CIF 构造，现用途映射将它们保持为同两种 basal termination 家族的 TRAIN 力标签候选；不能把其中任一点拆出作为 DEV。已有 48 原子体相及其位移数据属于同一 bulk-parent 家族。旧表面/界面候选有历史 MACE 训练来源、标签不完整，或派生自 135 原子界面；不构成合格的新 DEV 父族。既有 DFT 总表的父族/泄漏审计也指出，多份文件和不同 k 网格会重复表示同一几何，不能按文件数或网格数充当独立结构。

当前容器访问 `crystallography.net` 与 Materials Project 网页时均返回 `Tunnel connection failed: 403 Forbidden`。这只说明本容器当前无法取得外部结构，不代表这些数据库没有可用记录。

## 执行环境预检

| 项目 | 观测值 | 结果 |
|---|---:|---|
| CPU quota | 4 核（`400000 / 100000`） | 已读取 |
| cgroup 内存上限 | 34,359,738,368 字节（32 GiB） | 已读取 |
| `/workspace` 可用空间 | 16,643,678,208 字节 | 已读取 |
| Python `gpaw` 模块/命令 | 不存在 | 阻塞 |
| Python `ase` 模块 | 不存在 | 阻塞 |
| 当前 MPI/GPAW 作业 | 未发现 | 无作业冲突迹象 |

因此当前任务容器不能运行与既有标签一致的 GPAW/ASE 固定几何计算。本轮没有尝试安装软件；即使补齐运行依赖，也不能替代取得独立、可追溯的 DEV 结构源。

## 后续最小路径

1. 取得一个独立可追溯的 Ti₃SiC₂ 表面/结构父源（文件、原始来源 URL/文献、许可或来源说明及 SHA）；不要复用 COD 9009647 的派生几何作为“独立来源”。
2. 由该父源构建 DEV 参考结构和一个确定性 ±0.02 Å 局部响应结构，冻结原子映射、结构 SHA、输入 SHA 与同一方法锁；不改现有 TRAIN 候选。
3. 将两点送入具备 GPAW 26.7.0、ASE 3.29.0、匹配 PAW 数据且有持久结果归档的 E DFT 环境。当前 Codex 容器不具备这些软件。
4. TEST 结构身份、输入和标签继续由独立 evaluator 管理；E 不读取、不生成、不重标 TEST 数据。
5. 新 DEV 标签通过 SCF、几何映射和 SHA 校验后，才更新独立用途映射；本报告不改变任何数据角色。

若希望先做“不同晶面/表面家族”的内部 DEV，可从已有 COD 晶体源构造新的晶面 slab，但必须明确它与 TRAIN 共用同一个晶体学根源；它可测试未见表面取向，不能宣称为来源独立的 blind validation。本轮没有擅自采用这一降级定义。

## 证据文件

- `../tsc_backbone_dft_batch_20261017/PUBLIC_TASK_MANIFEST.json` SHA-256：`9997106fd1b2f8963453f9163ebbb98e4ab420e581d7c4488a4ca39d2c50b6ef`
- `../tsc_backbone_dft_batch_20261017/SURFACE_CANDIDATE_MANIFEST.json` SHA-256：`7a7a0cff495a7d2c87ec23dac813e21ee8f43a5ffc9c167f31afafcee77cdd30`
- `../nep_ace_dft_qualification_audit_20261015/DFT_MASTER_INVENTORY.csv` SHA-256：`48c8496ea2f88a97a865db41f6441fa018dfec3eaebbcd46cfad3cb4bcc36406`
- `../tsc_surface_train_label_qualification_20261020/SURFACE_LABEL_QA.json` SHA-256：`74464ed7f47f483d25ba29e87cbb37d450c5cb29a114c40bc5ccf284c2deb203`

上述为本地工作区扫描结果；本次未能对外部结构数据库完成在线检索。
