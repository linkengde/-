# E 临时工作区与 GitHub 恢复验收

采集时间：2026-10-11 01:58 UTC。范围仅为 Codex Cloud 临时工作区、非科学性 GitHub 恢复演练，以及 Ag64 暖机成本估算；未启动采样或其他计算。

## 临时工作区边界与容量

- 本任务工作目录为 `/workspace`，E 的独立仓库为 `/workspace/-`，分支为 `cloud-e-b1-dft`。本任务允许写入的路径边界是 `/workspace` 和 `/tmp`；`/workspace` 是可写 overlay，`/tmp` 是可写 tmpfs。各处均通过创建、fsync 并移除仅含测试文字的探针文件验证实际写入。没有写到边界外。
- `/workspace` 当前总容量 33,770,192,896 B，可用 28,798,566,400 B，可用 inode 2,039,418。`/workspace/scratch` 和 `/workspace/shared` 都落在同一 overlay，不是独立卷。`/tmp` 可用 17,585,635,328 B，但属于易失 tmpfs。未发现 `/mnt/e-work`、`/data` 或 `/scratch` 挂载。
- 依据用户提供的环境说明，本工作区按 Codex Cloud 的任务隔离临时环境处理；overlay 和 tmpfs 均不作为跨会话持久存储证据。无需申请普通云服务器挂载卷。
- 推荐的任务内候选目录为 `/workspace/.work/e/sampling`，当前未创建，等待 C 在 Codex 专用配置中签署路径后再建立。`DFT_WORK_ROOT` 不设置，DFT 继续 HOLD。

## GitHub 分支传输与恢复演练

- E 的独立分支 `cloud-e-b1-dft` 在本地和 GitHub 均指向 `6e517042eb5d75fb885f9c3cd886e342bf0547bd`；该提交是前一轮将 E 资源报告推送到 GitHub 的提交。本轮通过 GitHub Raw 从此分支下载 `STORAGE_EVIDENCE.md`、`RESOURCE_PREFLIGHT.json` 和 `SHA256SUMS.txt`，落到新的 `/tmp/e-github-restore.3F1Zaz/download`，`sha256sum -c` 两次均通过；将文件复制到独立 `recovered` 目录后再次校验通过。
- 这证明小型非科学性文件的独立分支上传、远端分支确认、下载、哈希校验和恢复流程可用。演练没有创建科学数据或改写历史文件。科研轨迹/检查点本身尚不存在，因此其实际归档恢复尚未测试。
- 环境有 Git LFS 3.6.1，但没有 LFS 文件，`git lfs env` 显示 `AccessUpload=none`、`AccessDownload=none`。本次只验证普通 Git 文件；大文件/LFS 配额和上传权限仍未验证，不将其视为已验收。

## Ag64 暖机归档与容量预算

暖机配置目前仍只是模板：1,000 步、0.5 fs、每 100 步写坐标帧、每 20 步写热力学行、最终生成检查点；每任务轨迹归档上限 536,870,912 B，检查点预留 33,554,432 B，最小空闲线 2,147,483,648 B。runner 会保存步 0 和步 1000 的重启 bundle，并在 bundle 中保存状态原子、压缩状态数组、轨迹、热日志、RNG 状态、配置/任务/来源哈希及内部 SHA 清单。

- 输出数量估算：每任务 11 帧（步 0、100…1000）和 51 条热日志（步 0、20…1000）。单帧仅 64 个原子的三维 float64 坐标基数为 1,536 B；再加 ASE 元数据、原子属性、热日志、步 0/步 1000 的重启副本和状态 JSON，预计未压缩任务目录低于约 1 MiB。E 两任务合计按 2 MiB 预留，远小于两任务各自 512 MiB 的提议上限。此为依据输出频率和数据形状的工程估算，未经实际轨迹测量。
- 建议归档只包含任务检查点、轨迹、热日志、来源/任务清单、性能记录和哈希清单；不把许可再分发状态未核实的 `Ag.eam.alloy` 参数文件打包。写入本地临时 `.partial` 归档并 fsync 后计算 SHA-256，再将归档和 sidecar 推到 E 分支；恢复时下载到新的空暂存目录，先校验外层 SHA 和 bundle 内部哈希、`task_id`、步数及 C 配置 SHA，再原子移到该任务的新恢复路径，绝不覆盖已有工作目录。完整采样如需每 1,000 步外传 checkpoint，归档频率、仓库增长和保留策略应由 C 签署。
- 当前 `/workspace` 余量在数值上足以容纳两任务约 2 MiB 暖机输出和 2 GiB 空闲门槛，但此余量只对当前任务会话有效，不能作为会话重启后的本地恢复保障。远端恢复可行性目前由小型非科学性文件演练证明；科学档案仍需在正式放行后逐项验证。

## 暖机时间估算与当前放行状态

- 没有 Ag64 采样性能实测，实际 1,000 步墙钟时间为 `UNVERIFIED`。仅作敏感性预算：若实际达到 0.1 s/步则约 100 s；若为 1.0 s/步则约 1,000 s。它们是速率情景，不是性能承诺。模板设置每任务暖机累计上限为 3,600 s；E 的两任务串行时累计最多按 7,200 s 规划，具体以 C 新签署配置为准。
- 读取的模板 `C_RELEASE_CONFIG.json` SHA-256 为 `57dc948d45daa3c51f7b3e5491a19e4f407546806531c3afae14a221b054ffc0`；`storage_release=false`、`warmup_release=false`、`continuation_release=false`。没有修改这些字段。

```text
OWNER=E
WORKSPACE_WRITE_TEST=PASS
GITHUB_SMALL_FILE_UPLOAD_DOWNLOAD_SHA256_RESTORE=PASS
SAMPLING_OUTPUT_ESTIMATE=<1 MiB/task; not measured
WARMUP_RUNTIME=UNVERIFIED; 100–1000 s scenario only; configured ceiling 3600 s/task
LOCAL_PERSISTENCE=NO (temporary overlay/tmpfs)
SAMPLING_READY=NO (await C's Codex-specific signed config and SHA)
DFT_READY=NO
```

本轮未安装软件、未删除既有数据、未创建采样/DFT 工作目录、未生成结构/轨迹、未运行 Ag64 模拟或 DFT。A/B 任务及其 checkpoint 未触及。
