# D 节点势函数与运行环境验收

核验日期：2026-10-11（Asia/Shanghai）。Owner 由当前分支 `cloud-d-b1-dft` 确认。本报告只包含文件身份、解析元数据、资源/身份预检与归档证据；未上传势文件、PLT 文件、论文或模拟数据。

## 验收状态

| 字段 | 结果 | 依据 |
|---|---|---|
| `LOCAL_POTENTIAL_PASS` | YES | 上传 ZIP 已在本机解压；文件字节数与冻结 SHA 完全匹配 |
| `ASE_EAM_PARSE_PASS` | YES | C 的 `inspect_eam.py` 使用现有 ASE 3.29.0 解析成功；仅读参数元数据 |
| `NODE_PREFLIGHT_PASS` | NO / HOLD | 所有资源、路径、机器指纹、本地势文件和既有远端归档检查通过；脚本唯一拒绝原因为 C 配置尚未签署 |
| `REMOTE_ARCHIVE_PASS` | YES | D 既有非科学合成归档已 push、fresh fetch、SHA 校验并恢复；本报告随后也会 fresh fetch 校验 |
| `BLOCKERS` | C release 未签署 | `CODEX_WARMUP_RELEASE_CONFIG.json` 的状态为 `C_HOLD_NOT_SIGNED_FOR_SCIENTIFIC_EXECUTION`，D `C_host_release=false`，`warmup_release=false` |

## 本机势文件身份与解析

- 输入 ZIP：`/workspace/attachments/66eb07de-ae8a-4b74-961f-079713369687/WMH2006_Ag_四个非PDF官方文件.zip`，SHA-256 `fd3670c0cf714c60b706c6b2d485ef6134b58b9fdb75aff86ca739ecf846d51c`。
- 解压目录：`/workspace/WMH2006_Ag-local-validation/`；没有替换或重新下载文件。
- 本机势文件：`/workspace/WMH2006_Ag-local-validation/Ag.eam.alloy`，780452 bytes，SHA-256 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`。ZIP 内四项 `SHA256SUMS.txt` 均通过。
- C 工具来源分支：`cloud-c-codex-ephemeral-warmup-20261011`，提交 `f174cf6adbe35101f993a613f725586e27de8bc1`。执行包 `SHA256SUMS.txt` 的 19 项全部通过，清单 SHA-256 为 `bc8f9f06e6e711e9b86f9d0814c76c6a1df8009795a74a398e657fc5ea18fd47`。
- `inspect_eam.py` SHA-256 `d3a437582dcd5bad7e361c6c81f64ac1191161170105b43133ffdd4065526a77`。首次调用系统 Python 时因该解释器没有 ASE 而未解析；改用既有 `/workspace/.venvs/gpaw-mpi/bin/python` 后通过，ASE 版本为 3.29.0。没有安装软件。
- 解析结果：`format=alloy`，元素 `Ag`，原子序数 47，cutoff `5.995011000293092 Å`，`nr=10000`、`nrho=10000`；解析器范围是参数表/插值元数据，不计算原子能量、力或轨迹。径向末端被工具标为 `C_TAIL_REVIEW_REQUIRED`（表格末端与 cutoff 相差一个 dr）；这不是能量/力验证，也没有改变当前 C HOLD。

## 节点预检

使用 C 配置 SHA-256 `207451d2d43279b9a1ac2d661bbca56b60641890781090e9ea14ecb55ebbd34a` 和 `node_preflight_codex.py`（SHA-256 `cbd10c7f969a0645cfa26446b3aa6a25dff9e4ed104492a17901a6e5dc6c46b5`），task `C-D-AGL-1250K-rho097-seed01`。原始工具 JSON 保存在同目录 `D_NODE_PREFLIGHT_RAW.json`，SHA-256 `dffd8945cdc7d61cb7c76811321965b4fe1b12ad6e7307e3c9b3e21305eb2cc7`。

当前快照：CPU quota 4.0 cores；cgroup memory limit `34,359,738,368` bytes；本进程 RSS `12,550,144` bytes；`/workspace` 可用 `28,794,232,832` bytes、可用 inode `2,038,456`。文件系统锚点为 `/workspace`。候选采样根 `/workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot` 不存在，archive staging `/workspace/.codex-ag64-archive/D` 与工作盘同文件系统；没有创建任何任务根。任务会话 fingerprint `6926836bb39914795d2a5304ad8f7326a8736f659508d110b1cc4e59584b1978` 与 C 当前配置所记录值相符。

除 `C_release_signed=false` 外，工具报告的 CPU、内存、空闲空间、inode、绝对路径、根目录无符号链接/不重叠/在 checkout 外、采样根缺失、本机势 SHA、会话 fingerprint 与既有合成归档测试全部为 `true`。因此资源证据已更新，但完整 `PREFLIGHT_STATUS` 仍为 `HOLD`，不是执行放行。

## GitHub 分支证据与限制

D 已有合格的非科学性往返测试：20,480-byte fixture SHA-256 `84165d5c3bdfaec96c774e6f2bd2429d6685261894b3c5aadfa2c21c7e4b737b`，上传提交 `5b599df646edd549f7a682974e170bce89fffd5e`；此前已从远端 fresh fetch、恢复到空目录并校验 archive 与全部成员 SHA。本轮按要求保留该证据，没有重复生成或覆盖 fixture。

本报告提交到 D 独立分支的新增提交 SHA 由 Git 提交结果及远端 fresh-fetch 校验记录。公开仓库只包含本报告和节点预检 JSON；`Ag.eam.alloy`、三个 PLT 文件、PDF 和其他势参数本体均未加入 Git。新增提交会使原 C HOLD 配置记录的 D 远端 HEAD 变旧；C 签署前需基于最新 D 分支头重新冻结配置 SHA。

本节点未启动 MD、DFT、AIMD、TTM 或训练；D 的其余采样、六任务续算和 DFT 继续 HOLD。
