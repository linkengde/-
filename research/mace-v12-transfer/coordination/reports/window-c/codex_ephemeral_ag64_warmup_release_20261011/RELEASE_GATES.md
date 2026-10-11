# 暖机放行门与停止条件

`CODEX_WARMUP_RELEASE_CONFIG.json` 是新模式配置模板，当前为 HOLD。它没有签署暖机：`config_status=C_HOLD_NOT_SIGNED_FOR_SCIENTIFIC_EXECUTION`、`release_decision=HOLD`、`storage_release=false`、`warmup_release=false`、三台 `C_host_release=false`，且每台 `potential_path=null`、没有本机势 SHA 收据。配置 SHA 仅表示字节完整性；启动命令必须另行由 C 发布精确 SHA，签署后任何修改都会使 SHA 不匹配并拒绝。

正式签署前，每个任务 owner 必须在其自己的 Codex 会话重跑 `node_preflight_codex.py`，上传报告并提供：当前 session fingerprint、根路径/同文件系统检查、实时可用字节和 inode、RSS/CPU quota/cgroup memory、准确的 `Ag.eam.alloy` 本机绝对路径和精确 SHA。势文件必须 SHA-256 等于 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`。配置不能下载替代势，也不能把 `persistence_verified=false` 改成 true。

远端归档先决条件：仓库/分支可达且可写；owner 独立分支的 live head 与配置相同；该 owner 已完成独立非科学性 fixture 上传、远端 fetch、哈希及恢复；可用空间至少满足 2 GiB 空闲线加 checkpoint 余量，inode 余量不少于 10,000；archive 单文件和总包均受 50 MiB 限制。E 目前 fixture 条件不通过。资源报告均为采集时快照，不当作下一次运行的实时测量。

恢复执行范围只允许逐任务暖机 1,000 步、0.5 fs，单任务墙钟最多 3,600 s，最大 RSS 8 GiB、空闲盘最低 2 GiB、checkpoint reserve 32 MiB；任务内部串行。参数、task IDs、seed、盒边、密度、温度路径、时间步和采样长度均来自原执行包，不允许本轮改动。只将 `C-D-AGL-1250K-rho097-seed01` 设为首个候选；本配置尚未授权执行。剩余五条暖机必须由 C 在首条归档和恢复审核通过后逐项更新签署，不能自动续算完整 120,000 步。

下列情况拒绝在创建采样根、初始化原子或计算力之前启动：配置未签署或 SHA 错误、模式不符、owner/task 错配、节点指纹变化、task 根重叠或位于 checkout 内、持久标志被改真、远端通道/fixture/分支头不合格、单文件或 archive cap 不合格、正确势文件/本机 SHA 收据缺失、资源预算不够、配置/科学参数与冻结基准不一致。运行中墙钟、RSS、磁盘、短程距离、非有限值、势函数表边界等任一 stop 条件触发即保存当前有效 checkpoint（若可行）并停止；不得静默调参。

本轮仍保持：A/B 不动；G/H 待命；其他 24 项、六个 GPAW 单点、DFT、完整采样续算、MACE 训练及 TTM 全部 HOLD。
