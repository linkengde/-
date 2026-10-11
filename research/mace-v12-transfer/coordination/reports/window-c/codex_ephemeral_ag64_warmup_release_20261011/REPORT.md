# Ag64 Codex 临时存储验收与暖机放行审查

**决定：`WARMUP_RELEASE = HOLD`。** Codex 专用存储模式、归档逻辑和拒绝执行测试已准备并核验；D/E/F 当前没有任何一个节点达到科学暖机放行条件。本次没有创建节点采样目录，没有运行结构生成、MD、DFT、LAMMPS、TTM 或 MACE 训练。

> 历史快照说明：以上状态是初始审计报告（配置 SHA `207451d2d43279b9a1ac2d661bbca56b60641890781090e9ea14ecb55ebbd34a`）。D/E/F 后续回执已补齐。当前仅 D 首项获得新配置的 1,000 步暖机条件放行；最新依据、签署及 SHA 见 `C_RELEASE_ATTESTATION_D_WARMUP_20261011.md`。本文件其余历史数据不重写。

首个审查候选仍为 D 的 `C-D-AGL-1250K-rho097-seed01`。这只是顺序指定，不是执行许可。不得自动启动同节点第二项、E/F 其他四项或余下 119,000 步。

## 节点证据

三个节点使用各自 Codex owner 分支及独立配置路径。分支文件表明归档 checkout/恢复的小型合成测试已在 D、F 通过；E 只验证了资源/存储报告文件的小型复制和哈希，没有独立 deterministic fixture。因此本配置将 E 的 `synthetic_fixture_verified` 保持为 false。D/E/F 均报告临时 overlay，不代表物理主机互斥，也不提供节点间硬件独占证据。

| Owner | 报告时 CPU 配额 | 报告时可用内存 | 报告时工作盘 | 可用 inode | 合成往返 | 节点本地精确势 SHA | 当前结论 |
|---|---:|---:|---:|---:|---|---|---|
| D | 4 CPU | 30,222,389,248 B | 28,796,854,272 B | 2,038,496 | PASS，20,480 B fixture | 未提供 | HOLD |
| E | 4 CPU | 30,511,386,624 B | 28,798,566,400 B | 2,039,418 | HOLD；只有报告文件 roundtrip | 未提供 | HOLD |
| F | 4 CPU | 30,780,076,032 B | 28,873,072,640 B | 2,036,714 | PASS，16,415 B fixture | 未提供 | HOLD |

数值是 owner 报告中的采样时读数，不是当前实时预检。暖机 runner 按原配置检查至少 2 GiB 空闲盘、32 MiB checkpoint reserve、10,000 个可用 inode、RSS 上限 8 GiB、单阶段 3,600 秒墙钟；这没有把 30 GB DFT 建议机械用作 Ag64 门槛。节点物理主机独占状态标为 `UNVERIFIED_NOT_CLAIMED`。提出的节点根路径在报告中不存在，均位于 `/workspace` checkout 外；它们须在签署时重新检查，不能把建议路径表述为已创建的持久目录。

详细证据提交、路径和文件 SHA 见 `OWNER_EVIDENCE_INDEX.csv`。提交号仅标识来源快照；不能替代下次运行前的新 owner 预检。

审查时从 GitHub 只读查询到 owner 归档分支头分别为 D `cloud-d-b1-dft`=`81fd10b429de56542135f2186de22291947e3189`、E `cloud-e-b1-dft`=`5b9ba5c3804df85d9910cc83ac30c4c29b58ac57`、F `cloud-f-b1-dft`=`8d4a721c6c4eb4588a61f994d325800d8f02b7d8`。这些头与本配置中的证据头相符；每个 owner 启动前仍必须重新查询，且普通 fast-forward 是唯一允许的写入方式。

## 放行阻塞

1. 三节点均没有可供 C 验证的本机 `Ag.eam.alloy` 文件路径及 SHA-256 回执。要求精确匹配 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`；程序不会下载或替换势函数。
2. E 缺独立的非科学合成 fixture 上传、fresh fetch、恢复及 SHA 测试。
3. 所有资源/节点指纹证据均需在启动前刷新；现有机器报告是临时容器的历史快照。
4. 远端分支写入权限、活动分支头、临时空间与归档暂存空间、根目录及同文件系统情况必须由各 owner 在同一会话重新核验。
5. 当前发布配置明确为 HOLD：不包含节点释放、存储释放或暖机释放。配置 SHA 是完整性标识，不是科学执行授权。

当前模板的 SHA-256 为 `207451d2d43279b9a1ac2d661bbca56b60641890781090e9ea14ecb55ebbd34a`。Codex runner SHA-256 为 `6b54d094c26bb9d4a5e4eb5997a31955eb8829a2975042e23fe10527cb943583`；完整的五程序 SHA 在配置的 `execution_code_sha256` 字段中。配置仍为 HOLD，不能用于启动。

## Codex 专用存储模式

新增模式 `CODEX_EPHEMERAL_WITH_VERIFIED_REMOTE_ARCHIVE` 明确保持 `persistence_verified=false`。独立配置、根路径、owner/task 范围、配置摘要、任务表摘要、节点会话指纹、远端分支及 head、归档通道、节点势文件 SHA 回执、资源约束和代码文件 SHA 在任何原子/计算器初始化之前检查。门槛不满足时退出；不会创建采样根。

暖机结束后，只有完整 checkpoint 原子写入并同步后，程序才创建任务归档。归档包括原子/坐标/速度/动量/质量/晶胞/PBC/ID，RNG 状态、温控与步数状态、restart、轨迹、热日志、来源和代码/任务/配置 SHA。每个文件和总归档均上限 50 MiB；超过即停止 Git 归档流程。程序先本地打包，推送 owner 独立分支，再以新 clone 取回，核验外层与内部文件 SHA、task/step/config 身份、恢复内容；上传恢复回执后再次 fresh fetch 验证。任何步骤失败都阻止后续 task 或续算。

公开仓库只接收用户授权的 checkpoint/轨迹和校验清单；不上传 `Ag.eam.alloy`、本 PDF 或许可未确认的论文全文。仓库公开可读已由 C 侧匿名 Raw HTTP 200 检查；运行时的 Git 写权限和分支 HEAD 仍须节点本地验证。

## 验证和范围

代码保留独立 `PERSISTENT_WORK_ROOT` 路径及其 `persistence_verified=true` 条件；Codex 临时模式拒绝把持久标志改成 true。纯合成测试覆盖：未签署配置、配置摘要错误、代码 SHA 错误、owner/task 错配、缺失归档通道、E fixture 未通过、势 SHA 错误/缺本机回执、持久模式兼容、临时模式持久标志约束、根路径创建前拒绝，以及本地 bare Git 归档 fresh fetch / SHA / 恢复。最后一次运行 `test_codex_release_gate.py` 的 17 项测试全通过；代码编译和 `--help` 检查通过，详见 `TEST_RESULTS.md`。

本目录不证明任何 Ag 构型已经生成或经过液态 QA。warmup 只可用来核实暖机执行、资源和归档流程；温度达到 1250 K 不能单独证明液态，也不能替代 RDF、Q6、MSD、独立 seed 和几何去重审核。

保持 A/B 任务与 checkpoint 不变；G/H 待命；余下 24 项、所有 GPAW DFT、MACE 训练和生产 TTM 继续 HOLD。
