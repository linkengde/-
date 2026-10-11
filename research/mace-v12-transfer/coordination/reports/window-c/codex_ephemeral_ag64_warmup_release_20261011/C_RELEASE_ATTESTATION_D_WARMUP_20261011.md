# C 放行签署：D 首个 Ag64 暖机

> 历史签署版本。D 的执行器在旧配置下第 0 步因未携带已审核的径向尾部策略而安全停止。此文件中的旧 SHA `4aee9d64859736f35459259625598f767d5f3e0d8b33bc71f0abb5baf31639e2` 已被取代，不要用于预检或启动。当前配置 SHA 和最新命令见 `C_RELEASE_ATTESTATION_D_WARMUP_TAIL_REVIEW_20261011.md`。

日期：2026-10-11 UTC
状态：`APPROVE_WARMUP_ONLY`，范围严格限定为一个 task。

## 已核验依据

- D 当前 owner 分支头：`cloud-d-b1-dft` / `2c8b6f41ec8a7605ddc06f7218014c55957166ac`。报告 [`D_NODE_POTENTIAL_RUNTIME_ACCEPTANCE_20261011.md`](https://github.com/linkengde/-/blob/2c8b6f41ec8a7605ddc06f7218014c55957166ac/research/mace-v12-transfer/coordination/reports/window-c/five_cloud_preflight_and_liquid_pilot_readiness_20261011/owner_reports/D_NODE_POTENTIAL_RUNTIME_ACCEPTANCE_20261011/D_NODE_POTENTIAL_RUNTIME_ACCEPTANCE_20261011.md) 和原始 preflight JSON 已从该分支读取。报告 SHA-256 为 `6ce619eea8da27e3aff97a93a3893d09187cc50eb2a41410466a051072f8d2f7`；原始 JSON SHA-256 为 `dffd8945cdc7d61cb7c76811321965b4fe1b12ad6e7307e3c9b3e21305eb2cc7`。
- D 的 780452-byte 本地 `Ag.eam.alloy` SHA-256 与冻结值 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447` 一致；ASE 3.29.0 metadata-only 解析通过。势文件本体未上传。
- D 最近 preflight 的路径、身份、CPU、内存、磁盘、inode、本机势和既有远端归档子检查均通过；总状态仅因当时 C 配置未签署而为 HOLD。D 报告快照显示 4 CPU quota、32 GiB cgroup memory limit、28,794,232,832 free bytes、2,038,456 free inodes；路径根尚未创建，持久标志保持 false。
- E 最新 owner 分支和专用合成往返回执已由 C 核验：分支头 `39734568216e5a3d9853d7d052d28f9799172598`，报告清单 7 项均通过；确切势 SHA/ASE metadata parse PASS、20,480-byte fixture 的 fresh fetch/restore/byte comparison/SHA PASS。E 不获本次执行放行。
- F 最新 owner 分支已核验：分支头 `2ddb39a899decd6b3cc8d6c002d48291c9a5e632`，报告清单 6 项均通过；确切势 SHA/ASE metadata parse PASS，当前两项资源预检子检查通过，既有 synthetic roundtrip PASS。F 不获本次执行放行。
- 候选密度来源及 Ag64 盒尺寸已核验。cutoff 大于盒半边长，但现有来源审查确认本执行器固定使用 ASE `PrimitiveNeighborList` 的整数周期 offset 并保留同索引的不同周期像。该支持只解决邻居枚举，不能证明 64 原子盒的熔点、扩散或凝固物性。
- Codex 临时归档通道限制仍固定：每个文件及归档不超过 50 MiB，运行后必须 push、fresh fetch、恢复并核对内外层 SHA；任何备份失败均停止，不得续算。不得上传势文件、PLT、PDF或论文全文。

## 签署配置

现行 `CODEX_WARMUP_RELEASE_CONFIG.json` SHA-256：

`4aee9d64859736f35459259625598f767d5f3e0d8b33bc71f0abb5baf31639e2`

配置保持冻结的势、任务矩阵、密度/盒子、seed、温度、步长和采样参数；`persistence_verified=false`。它只将 D 的 `C_host_release` 设为 true，并仅将 `C-D-AGL-1250K-rho097-seed01` 列入 `warmup_task_ids`。`continuation_release=false`，没有自动授权该任务剩余 119000 步；E/F host release、D 第二个 task、E/F 其余 tasks、DFT、AIMD、TTM及MACE训练均保持 HOLD。

这是 C 在仓库提交中发布的授权声明，配置 SHA 和 Git commit 用于固定可复核字节；这不是公钥数字签名。

## D 节点执行前最后一道本机门

D 必须从 C 分支取到此签署配置及匹配的执行包，在自己的当前 Codex 会话中先运行只读预检，并将 task 根保持未创建：

```bash
python node_preflight_codex.py \
  --config CODEX_WARMUP_RELEASE_CONFIG.json \
  --config-sha256 4aee9d64859736f35459259625598f767d5f3e0d8b33bc71f0abb5baf31639e2 \
  --owner D \
  --task-id C-D-AGL-1250K-rho097-seed01 \
  --potential-path /workspace/WMH2006_Ag-local-validation/Ag.eam.alloy
```

只有此命令返回 `PREFLIGHT_STATUS=PASS` 且运行时没有新变化时，才按原有用户授权运行：

```bash
python run_sampling_codex.py \
  --config CODEX_WARMUP_RELEASE_CONFIG.json \
  --config-sha256 4aee9d64859736f35459259625598f767d5f3e0d8b33bc71f0abb5baf31639e2 \
  --owner D \
  --task-id C-D-AGL-1250K-rho097-seed01 \
  --stage warmup
```

执行器会再次检查配置 SHA、代码 SHA、D owner/task、远端分支头、资源、势 SHA及存储门。任何一项非 PASS 即停止；不得改参数绕过。不得在运行前先向 D 分支提交别的文件，因为这会使配置冻结的远端 HEAD 失效。

成功时程序最多执行 1000 步 × 0.5 fs，自动执行规定的 checkpoint/归档/远端恢复校验，并将状态标为待 C 审核；不自动进入完整采样。D 完成后应报告新分支提交、恢复 SHA、实际资源和墙钟；C 审核通过前不启动 D 的第二个 task。
