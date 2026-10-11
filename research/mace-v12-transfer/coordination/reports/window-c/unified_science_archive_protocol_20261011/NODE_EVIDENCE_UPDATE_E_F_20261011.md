# E/F 节点非科学性验收证据更新

日期：2026-10-11
用途：更新 C 的统一节点证据索引；不签署科学执行配置。

## 结论

E 和 F 均已在各自 owner 分支提交当前 Ag 势文件身份、ASE 元数据解析、Codex 资源预检及 GitHub 远端归档证据。C 从各自远端分支读取原始报告并按报告内 `SHA256SUMS.txt` 独立校验，所有列出的文件哈希通过。两个节点仍未获科学运行放行：C 配置 `207451d2d43279b9a1ac2d661bbca56b60641890781090e9ea14ecb55ebbd34a` 仍为 HOLD，预检也如实返回 HOLD。

## E 节点

- owner 分支：`cloud-e-b1-dft`；核验分支头：`39734568216e5a3d9853d7d052d28f9799172598`。
- 验收目录：`research/mace-v12-transfer/coordination/reports/window-c/liquid_ag_six_pilot_authorized_dispatch_20261011/owner_reports/E/codex_node_acceptance_20261011/`。
- `ACCEPTANCE_REPORT.md` SHA-256：`cf3dfe60a966963705e34c2eacee949b1c436e4ab6082d18a14e25d70e9c29c6`；该目录七项清单校验均通过。
- node preflight JSON SHA-256：`bc58ce061313e4319038cd72b6fa1607f90b8b01f23092dd7925096c8f9b382a`。报告中记录 E 本机 `Ag.eam.alloy` 为 780452 bytes、SHA-256 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`；ASE 3.29.0 元数据解析通过，没有计算势能、力或运行 MD。
- 远端合成往返：20,480-byte 独立合成 fixture，SHA-256 `c3c582d73df54e9e9118efec0b96ca1c8faebba7b363ce3ebb007796cb250cf6`；全新 Git 仓库 fetch、恢复、字节比较及 SHA 校验 PASS。fixture 上传提交为 `649b7f4ef673ac618cf962f804f8081ff993ffc5`；当前 E 分支头为后续验收报告提交 `397345...`。
- 节点预检结果为 HOLD，原因包括旧 C 配置仍未签署、当时配置中的 synthetic-roundtrip 状态仍为 HOLD。CPU quota 4、cgroup memory limit 34,359,738,368 bytes、工作区可用 28,795,498,496 bytes；Codex 工作区持久性仍为 false，物理主机独占性未声称。
- cutoff 几何核查：E 解析摘要显示 cutoff `5.995011000293092 Å` 大于登记 Ag64 盒半边长。C 已归档的 `SAMPLER_SOURCE_REVIEW.md` 明确核对了实际采用的 ASE `PrimitiveNeighborList` 整数周期 offset：它会保留同一原子在不同周期像中的邻居，因此支持 cutoff 大于半盒长的多周期像；执行器不能将邻居按原子索引去重，也不能改成最小像近似。该问题不是额外放行阻塞。此核查不验证 64 原子有限尺寸下的熔点、扩散或凝固物性；这些仍需后续独立验证。

## F 节点

- owner 分支：`cloud-f-b1-dft`；核验分支头：`2ddb39a899decd6b3cc8d6c002d48291c9a5e632`。
- 验收目录：`research/mace-v12-transfer/coordination/reports/window-c/five_cloud_owner_evidence_20261011/F/codex_cloud_temp_storage_acceptance_20261011/ag_potential_runtime_preflight_20261011/`。
- `ENVIRONMENT_ACCEPTANCE.md` SHA-256：`d777b14fd440917698521de63a0cdb0b27131b0dbf178936407ba0fd6773c5ee`；该目录六项清单校验均通过。
- 势文件身份为 780452 bytes、SHA-256 精确匹配冻结值；ASE 3.29.0 metadata-only parse 通过，没有进行能量/力计算。原有 one-bin radial-tail 策略记录为已审查；此次仍需独立解决上述 cutoff/半盒长的周期性适用问题。
- 两个 F task 的 node preflight JSON 分别为 `C-F-AGL-1250K-rho100-seed02_node_preflight.json` 和 `C-F-AGL-1250K-rho103-seed01_node_preflight.json`，SHA-256 分别为 `3edd2a000c4708f5fb942e9f594c292ebf86ae801ae41d7b06115c28615f30d7` 与 `58099c5694a1879c4915a3268f5c01ed40ed338ff16e1110fee2488345627833`。资源、节点指纹、路径、势文件和归档子检查均通过；整体 HOLD 唯一原因是 C/owner release 尚未签署。
- 报告保留既有合成恢复测试 PASS，并记录本次报告上传后从远端 fresh fetch 的 SHA 校验。工作区约 28.87 GB 可用空间、约 2.04M 可用 inode；工作区为临时空间，持久性未验证，不能将其表述为持久盘。

## 统一状态与下一步门槛

| 节点 | 势文件身份/解析 | 远端归档 | 节点预检 | 科学任务 |
|---|---|---|---|---|
| E | PASS / PASS | PASS | HOLD（原 preflight 使用旧配置；新配置签署后需刷新） | HOLD |
| F | PASS / PASS | PASS | HOLD（C 配置未签署） | HOLD |

这两份回执没有授权生成结构或运行 MD、DFT、MACE 训练。C 已归档的采样器来源审查解决了 cutoff 多周期像处理；这不改变液态物性验证要求。C 配置只允许 D 的首个暖机 task ID，E/F 仍不能因环境验收通过而启动。A/B 不变，G/H 待命，其余 24 个 DFT 计划项保持 HOLD。
