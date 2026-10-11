# Codex 临时工作区归档协议

本目录中的新模式为 `CODEX_EPHEMERAL_WITH_VERIFIED_REMOTE_ARCHIVE`。它承认 `/workspace` 是任务临时空间，不把临时盘写成持久盘。配置字段 `persistence_verified` 始终为 `false`。旧 `PERSISTENT_WORK_ROOT` 模式和它的 `persistence_verified=true` 检查留在原 runner 路径中。

每个任务使用独立根：D `/workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot`，E `/workspace/.work/e/sampling`，F `/workspace/.work-f/sampling`。归档暂存分别放在 `/workspace/.codex-ag64-archive/{D,E,F}`。这些目录都在 `/workspace` 可写范围内、位于 checkout 外，彼此不重叠；它们表示路径隔离，不证明各 Codex 任务使用物理独占主机。当前它们只是待签署路径，报告显示尚未创建。运行前必须在 owner 节点重新确认实际路径、父目录挂载、空闲字节、inode、RSS、CPU 配额、内存限制、同一文件系统和根目录不存在；runner 在创建采样目录前再次核实采样/归档空间、inode 及设备号。DFT 工作根与这组采样根分开；DFT 仍 HOLD。

临时模式配置钉住 `run_sampling_codex.py`、`inspect_eam.py`、`codex_release_gate.py`、`archive_checkpoint.py` 和 `node_preflight_codex.py` 的 SHA-256。runner 在任何原子初始化之前比对实际程序字节；代码清单错误或不全即拒绝。原始持久模式脚本仍以其原始 SHA 保留，且持久模式继续要求 `persistence_verified=true`；新临时模式要求它保持 false。

暖机完成后，runner 只在最后完整 checkpoint 原子写入并 `fsync` 后打包。包必须含 `sampler_state.json`、最终 restart bundle（原子状态、坐标、速度、动量、质量、晶胞、PBC、原子 ID）、PCG64 各 RNG 状态、温控参数和状态、完整 `trajectory.traj`、`thermal.csv`、`source_manifest.json`、任务表 SHA、配置 SHA、代码 SHA、性能日志、budget 记录及内部 SHA 清单。restart bundle 的文件 SHA 必须与状态指针一致。

本地临时归档先写 `.partial`、`fsync` 后原子改名。每个文件和压缩档案都不得大于 50 MiB（52,428,800 bytes）；任何单文件或档案超过上限即停止 Git 归档路径。GitHub 分支采用 owner 独立分支和新的任务专属目录，使用普通 fast-forward 提交，绝不 force push；上传前远端分支必须仍等于签署配置中的 expected head。只上传任务 checkpoint 与上述清单，不上传 `Ag.eam.alloy`、PDF、论文全文或其他未确认再分发许可的文件。

上传后 runner 新建临时 Git clone，从远端该提交重新获取压缩档案和清单，先核对外层 SHA，再检查安全 tar 路径、成员文件集合、每个内部 SHA、task ID、step=1000、配置 SHA、任务表 SHA和restart checkpoint 状态。通过后生成 readback receipt，再以另一个普通提交上传 receipt，最后再次从远端 fresh clone 拉回档案和 receipt并校验。两次提交中任一步、上传、读取、恢复、哈希或身份核验失败，runner 返回 HOLD；不允许续算或继续下一任务。临时 staging 在结束时清理，owner task 根和 GitHub 远端保留任务输出及恢复证据。

GitHub 仓库匿名 Raw 读取返回 HTTP 200，仓库为公开可读。用户指令明确要求将这些 checkpoint/轨迹经独立分支做可验证远端备份，因此本方案仅上传用户授权的科学输出；不上传模型文件或许可不明资料。分支写权限、归档保留和远端实时分支头仍须由每节点在运行前检查。超过 50 MiB 的文件/归档、分支冲突、上传或回读失败均停止，改用另行批准的存储通道前不得继续。

此协议的 owner 侧小型合成演练目前 D 和 F 通过，E 仍缺独立合成载荷演练。C 的本地 bare-Git 集成测试只验证归档程序实现，不替代 D/E/F 各自的 Codex 环境验收，更不证明生产输出可长期保存。
