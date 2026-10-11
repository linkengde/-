# 暖机放行门与停止条件

当前 `CODEX_WARMUP_RELEASE_CONFIG.json` SHA-256 为 `4aee9d64859736f35459259625598f767d5f3e0d8b33bc71f0abb5baf31639e2`。C 只签署 D 的 `C-D-AGL-1250K-rho097-seed01`、`warmup` 阶段，最多 1,000 步、0.5 fs；D host release 为 true，E/F host release 为 false。`continuation_release=false`，其余五个暖机 task、剩余 119000 步、全部 DFT、MACE 训练及 TTM 不在本次放行范围。签署前原始 HOLD 配置保存在 `CODEX_WARMUP_RELEASE_CONFIG_PRE_RELEASE_HOLD.json`。精确配置 SHA 是启动必要条件，签署后改动任何配置字节都会令 SHA 失效。

D 的最新 owner 回执显示本机势、资源、路径、会话指纹和归档子检查通过；整体 HOLD 原因是当时 C 配置未签署。D 必须在自己的当前 Codex 会话以新 SHA 重跑只读 `node_preflight_codex.py`，仅当结果为 PASS 才启动 runner；runner 在生成坐标前还会独立重检相同条件。势文件 SHA-256 必须等于 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`。配置不能下载替代势，也不能把 `persistence_verified=false` 改成 true。E/F 最新回执也已核验，但其 host release 仍为 false，不可启动。

远端归档先决条件：仓库/分支可达且可写；owner 独立分支的 live head 与配置相同；该 owner 已完成独立非科学性 fixture 上传、远端 fetch、哈希及恢复；可用空间至少满足 2 GiB 空闲线加 checkpoint 余量，inode 余量不少于 10,000；archive 单文件和总包均受 50 MiB 限制。D/E/F fixture 当前均 PASS。资源报告都是采集时快照，D 启动前必须按签署配置 SHA 在当前会话重测。

本次执行范围只允许 D 首项暖机最多 1,000 步、0.5 fs，单任务墙钟最多 3,600 s，最大 RSS 8 GiB、空闲盘最低 2 GiB、checkpoint reserve 32 MiB；任务内部串行。参数、task ID、seed、盒边、密度、温度路径和步数保持原执行包冻结值。D 其余任务、E/F 全部任务和另外 119,000 步不得自动续算。后续五项暖机必须由 C 在审核首项实际归档和恢复结果后逐项另签。

下列情况拒绝在创建采样根、初始化原子或计算力之前启动：配置未签署或 SHA 错误、模式不符、owner/task 错配、节点指纹变化、task 根重叠或位于 checkout 内、持久标志被改真、远端通道/fixture/分支头不合格、单文件或 archive cap 不合格、正确势文件/本机 SHA 收据缺失、资源预算不够、配置/科学参数与冻结基准不一致。运行中墙钟、RSS、磁盘、短程距离、非有限值、势函数表边界等任一 stop 条件触发即保存当前有效 checkpoint（若可行）并停止；不得静默调参。

本轮仍保持：A/B 不动；G/H 待命；其他 24 项、六个 GPAW 单点、DFT、完整采样续算、MACE 训练及 TTM 全部 HOLD。
