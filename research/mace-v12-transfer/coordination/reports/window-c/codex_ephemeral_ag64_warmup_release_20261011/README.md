# C：Ag64 Codex 暖机存储验收

本目录记录 D/E/F 六条 Ag64 候选任务的临时工作区归档门审查。执行模式承认 `/workspace` 为临时空间，`persistence_verified` 保持 `false`；旧持久工作根模式及其检查另行保留。

当前 `CODEX_WARMUP_RELEASE_CONFIG.json` 已由 C 补充签署，但只条件放行 D 的首项 `C-D-AGL-1250K-rho097-seed01` 暖机，最多 1,000 步、0.5 fs；D 必须先用最新签署 SHA `2e122b258a48fe54bf603709b5505a301c1257f398c5625eecc2ca94998c0052` 完成本机只读 preflight 且获得 PASS。E/F 的 `C_host_release=false`，D 第二任务和其余四条均 HOLD；continuation、完整采样、DFT、训练和 TTM 均未放行。已签但缺少径向尾部策略字段的旧配置字节保存在 `CODEX_WARMUP_RELEASE_CONFIG_PRE_TAIL_POLICY_20261011.json`；最初 HOLD 配置保存在 `CODEX_WARMUP_RELEASE_CONFIG_PRE_RELEASE_HOLD.json`。

D/E/F 均已提供本机势 SHA、ASE metadata-only 解析和当前非科学性归档证据，C 已在 `OWNER_EVIDENCE_INDEX.csv` 登记相应 owner 分支及文件 SHA。Codex 临时目录仍保持 `persistence_verified=false`，也不声称底层物理主机互斥。

执行顺序：先查看 `D_WARMUP_RELEASE_SUMMARY_20261011.md`、`C_RELEASE_ATTESTATION_D_WARMUP_TAIL_REVIEW_20261011.md`、`NODE_RELEASE_MATRIX.csv`、`OWNER_EVIDENCE_INDEX.csv` 和 `RELEASE_GATES.md`。D 节点先用最新 attestation 中的 SHA 运行 `node_preflight_codex.py`，仅当输出 PASS 才能启动指定 1,000 步暖机。runner 会再核验配置/代码/节点/势/资源/分支头/归档门；任一失败即拒绝。径向尾部审查仅覆盖原样势文件 SHA、ASE 3.29.0 和单 `dr` 插值范围。

`test_codex_release_gate.py` 包含本地合成数据与本地 bare Git 远端测试；它验证程序拒绝条件和归档代码的读回恢复流程，不替代 D/E/F 节点验收，也不产生材料模拟数据。`AUTHORIZED_SIX_TASKS.csv` 与持久模式配置/runner的基准副本只用于证明六任务映射及原持久模式保持不变。

本目录不含 Ag.eam.alloy 文件或论文 PDF。远端归档协议禁止上传许可未确认的势文件、论文全文及其他未获再分发许可的资料。

当前运行授权、owner 最新证据和执行前命令见：

- `D_WARMUP_RELEASE_SUMMARY_20261011.md`
- `C_RELEASE_ATTESTATION_D_WARMUP_TAIL_REVIEW_20261011.md` (latest; use this SHA)
- `NODE_RELEASE_MATRIX.csv` 和 `OWNER_EVIDENCE_INDEX.csv`
- `CODEX_WARMUP_RELEASE_CONFIG_PRE_RELEASE_HOLD.json` 保存签署前原始配置字节
