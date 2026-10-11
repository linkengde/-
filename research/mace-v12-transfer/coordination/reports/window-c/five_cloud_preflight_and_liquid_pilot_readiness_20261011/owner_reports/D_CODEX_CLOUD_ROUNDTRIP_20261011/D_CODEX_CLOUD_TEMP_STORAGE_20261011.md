# D Codex Cloud 临时存储替代验收

核验时间：2026-10-11 01:59 UTC。节点身份由独立分支 `cloud-d-b1-dft` 确认。本报告只验证 Codex Cloud 临时工作区及 GitHub 小文件往返恢复；没有启动 Ag64 采样、DFT、MD、AIMD、TTM 或训练。

## 结论

| 项目 | 结果 | 说明 |
|---|---|---|
| `/workspace` 可写 | YES | 工作区、仓库与 `.git` 的写入探针通过 |
| 任务目录边界 | VERIFIED | 当前仓库为 `/workspace/-`，位于授权可写的 `/workspace` 下；本次只改 D 分支 |
| 跨会话本地持久性 | NO / UNVERIFIED | `/workspace` 是 Kata overlay；没有平台跨任务保留或恢复证明。`/tmp` 是 tmpfs |
| GitHub 小文件上传、下载、SHA 与恢复 | PASS | 20,480-byte 非科学性 tar 在 D 分支往返并恢复到新目录，归档和成员哈希一致 |
| `SAMPLING_STORAGE_READY` | NO / HOLD | 远端小文件往返已证实，但尚无 C 对此替代存储方案的签署；未创建采样根 |
| `DFT_STORAGE_READY` | NO / HOLD | 未验证持久 DFT 根，也没有冻结输入后的 GPAW checkpoint 实测与保留预算 |
| `RESOURCE_READY` | NO / HOLD | overlay 容量可用，但本地跨会话持久性及 C 的存储放行仍缺失 |

## 工作区与容量

`findmnt -T /workspace` 显示 `overlay[/workspace]`、`overlay`、`rw`；`stat` 显示 `/workspace`、`/workspace/-` 和 `/workspace/-/.git` 位于同一设备。对 `/workspace`、当前 checkout、Git 元数据和 `/tmp` 的可写探针均通过。当前仓库是独立 D checkout，分支 `cloud-d-b1-dft`；没有访问或改动 A/B 工作区。

2026-10-11 01:59 UTC 的文件系统读数：`/workspace` 总计 `33,770,192,896` bytes，可用 `28,796,854,272` bytes；inode 可用 `2,038,496`。这些是 overlay 的可见容量，不是独立卷配额。没有创建 `/workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot` 或 DFT 根，也没有删、移动或压缩历史数据。

可在 C 后续签署后考虑的采样根为 `/workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot`。它仍在临时 overlay 内，现阶段仅是路径提案；创建该目录本身也不能证明持久性。DFT 仍使用既有提案 `/workspace/dft-work/cloud-d/ag-ti-si-c-b1`，但保持未创建、未设置和 HOLD。

## GitHub 往返与恢复测试

测试只包含通用说明文本和确定性字节序列，不含原子结构、科学参数、轨迹或 checkpoint。归档提交在 D 独立分支 `cloud-d-b1-dft`：

- 测试归档：`D_NONSCIENCE_ROUNDTRIP_FIXTURE.tar`，20,480 bytes，SHA-256 `84165d5c3bdfaec96c774e6f2bd2429d6685261894b3c5aadfa2c21c7e4b737b`。
- 上传提交：`5b599df646edd549f7a682974e170bce89fffd5e`；推送未使用 force。远端分支 fetch 后指向同一提交。
- 使用远端分支的 Git blob 下载归档，下载件 SHA 与上传前一致；解压到全新 `/tmp/D_NONSCIENCE_ROUNDTRIP_RESTORE_20261011` 后，`README.txt`、`SHA256SUMS.txt`、`payload.bin` 三个成员均通过 SHA-256 校验。
- 该结果证明当前认证与分支可传输、恢复 20 KiB 级普通文件；不证明临时 overlay 会跨会话保留，也不证明 GitHub 对 512 MiB 采样归档的传输、容量或长期仓库增长可行。

### 建议的 Ag64 状态归档与恢复流程

在收到 C 签署的 Codex 专用暖机配置及精确 SHA 前，不创建任务根、不打包科学数据、不上传 checkpoint。获准后，针对每个 D task 使用各自独立目录；只在 sampler 已完成并 fsync 的 checkpoint 边界归档，保留原始输出。归档应包含当前 `state.json`、其引用的完整 restart bundle（`state_atoms.traj`、`state_arrays.npz`、`trajectory.traj`、`thermal.csv`）、`source_manifest.json`、对应 task 的日志，以及绑定的配置/输入提交标识。归档外记录 SHA-256 清单和归档 SHA；每次生成新、不可覆盖的文件名，推送到 `cloud-d-b1-dft` 的独立路径，不改写 `main` 或其他 owner 分支。

恢复时，从记录的确切分支提交读取归档到新建空目录，先校验归档 SHA，再检查成员路径并校验完整 SHA 清单，确认 task ID、checkpoint step、输入 manifest SHA 和 C 配置 SHA 均匹配。只有 C 发布的 `resume_checkpoint_sha256` 与恢复状态一致后，才可按其后续放行执行续算。此流程仍需 C 接受“GitHub 分支作为远端备份”的存储角色；当前测试未上传任何科学数据。

C 的每任务轨迹 archive cap 为 512 MiB，高于普通 GitHub 单文件约 100 MiB 的限制；20 KiB 测试不能外推到大归档。若任一归档预计超过该传输限制，或累计 Git 仓库增长未获 C 接受，应保持 HOLD 并由 C 指定合适的数据归档通道。完整 120,000 步轨迹和 DFT checkpoint 不因此获得放行。

## 1000 步暖机预算估算

以下是静态代码审查估算，不是已运行或测量的结果。固定执行包 `2959fd3e2396a7319bd1f5c7f4fec27ce442abe0` 的 `run_sampling.py` SHA-256 为 `b78b734798c9c075826e6611e3d2f7f156f0b05ff7c87159812ac76ed0e918a8`；固定 `C_RELEASE_CONFIG.json` SHA-256 为 `57dc948d45daa3c51f7b3e5491a19e4f407546806531c3afae14a221b054ffc0`。配置要求 Ag64、1000 warmup steps、每 100 步写一帧、每 20 步写热日志；初始化和暖机结束各保存一次 restart bundle。按 11 帧、51 行热日志，以及 bundle 内的轨迹/日志副本估算，单 task 输出约 `0.1–0.3 MiB`；暂按 `1 MiB/task` 预留，等待获准后用实际 `archive_bytes` 取代估算。该估值不构成 512 MiB archive cap 或 2 GiB 开跑前余量的验收。

固定配置给 warmup 阶段累计上限 `3,600 s/task`。由于没有获准的基准运行，本节点无法给出可靠的实测时长或更窄的预测区间；可用于资源排期的值为最多 1 小时/task，实际秒数和 steps/s 保持 `UNVERIFIED`，应在 C 签署并授权后由该 1000 步阶段记录。该时间上限本身不是开跑授权。

## 放行与保留项

固定配置中 `storage_release=false`、`warmup_release=false`、D 的 `C_host_release=false` 保持不变。既有资源报告与原始 preflight JSON 均保留，哈希分别为 `43815a97c18f865a433d234bd06fb12cefc3f20f1e9f59b20643a3482b3c393e` 与 `8be2a47f6dd17b94a2252cd72f9e379491ca46131a26f0a6db2b5944c13e73e7`。既有 `RUNTIME_PASS=YES` 结论保留，本轮未重跑软件或 MPI 检查。

当前 `PERSISTENCE_VERIFIED=NO`、`SAMPLING_STORAGE_READY=NO`、`DFT_STORAGE_READY=NO`、`RESOURCE_READY=NO`。等待 C 审核临时存储替代证据，并发布新的 Codex 专用暖机放行配置及 SHA；收到前暖机、完整采样续算、所有 DFT 均保持 HOLD。A/B 不受影响。
