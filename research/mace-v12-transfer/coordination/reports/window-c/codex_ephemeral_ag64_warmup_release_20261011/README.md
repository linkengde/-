# C：Ag64 Codex 暖机存储验收

本目录记录 D/E/F 六条 Ag64 候选任务的临时工作区归档门审查。执行模式承认 `/workspace` 为临时空间，`persistence_verified` 保持 `false`；旧持久工作根模式及其检查另行保留。

`CODEX_WARMUP_RELEASE_CONFIG.json` 是带精确代码 SHA 的 HOLD 配置，不是暖机批准。D/F 的小型非科学 fixture 已由各节点报告完成上传、fresh fetch、恢复及 SHA 核对；E 仅有报告文件的来回读取记录，没有独立 deterministic fixture，因此 E 的 fixture gate 仍 HOLD。D/E/F 均缺少可由 C 核对的节点本地 Ag.eam.alloy 精确 SHA 收据，节点资源数据也只是历史快照。没有节点因此获准启动。

执行顺序：先查看 `REPORT.md`、`NODE_RELEASE_MATRIX.csv`、`OWNER_EVIDENCE_INDEX.csv` 和 `RELEASE_GATES.md`；在各 owner 会话运行只读 `node_preflight_codex.py`；只有 C 更新配置、确认节点回执、冻结完整配置 SHA 并明确签署后，`run_sampling_codex.py` 才可能接受暖机。初次候选只限 D 的 `C-D-AGL-1250K-rho097-seed01`，1,000 步、0.5 fs。此配置目前拒绝所有科学执行。

`test_codex_release_gate.py` 包含本地合成数据与本地 bare Git 远端测试；它验证程序拒绝条件和归档代码的读回恢复流程，不替代 D/E/F 节点验收，也不产生材料模拟数据。`AUTHORIZED_SIX_TASKS.csv` 与持久模式配置/runner的基准副本只用于证明六任务映射及原持久模式保持不变。

本目录不含 Ag.eam.alloy 文件或论文 PDF。远端归档协议禁止上传许可未确认的势文件、论文全文及其他未获再分发许可的资料。
