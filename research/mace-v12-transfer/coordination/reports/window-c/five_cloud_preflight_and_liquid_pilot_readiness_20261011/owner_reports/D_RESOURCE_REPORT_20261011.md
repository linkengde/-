# D 云端持久工作区与资源补证

测量时间：2026-10-11 09:50:19 Asia/Shanghai（2026-10-11 01:50:19 UTC）。Owner 由当前独立分支 `cloud-d-b1-dft` 识别为 D。报告提交以本分支当时的 `f9544cd00046e7fd4be75fec4d6fe13e50e8cdab` 为父提交；它只是分支历史定位，不是本资源报告 SHA。

## 验收状态

| 字段 | 状态 | 依据 |
|---|---|---|
| `OWNER` | D | 本机分支为 `cloud-d-b1-dft` |
| `RUNTIME_PASS` | YES | 脚本语法、setup 路径、版本、PAW SHA 及四 rank GPAW MPI 初始化通过 |
| `PERSISTENCE_VERIFIED` | NO | 唯一可见工作盘 `/workspace` 为 overlay；没有 provider 持久卷策略或恢复证据 |
| `SAMPLING_STORAGE_READY` | NO | 可用字节足够满足 C 提案的小型采样预算，但没有合格的持久化独立根 |
| `DFT_STORAGE_READY` | NO | 持久根缺失，且没有冻结输入后的实际 checkpoint 大小和保留预算 |
| `RESOURCE_READY` | NO | 持久性、独立性及有效工作根尚未核实 |

## 固定提交与执行包

- 输入固定提交：`2959fd3e2396a7319bd1f5c7f4fec27ce442abe0`。
- 已读取该提交中的 `liquid_ag_six_pilot_authorized_dispatch_20261011/START_D.md`；文件 SHA-256 `d29c4a52f75c569aad4b02eba84925541e1d203e3f0fa2fa6b59c9a297d75619`。
- 对该执行包运行 `sha256sum -c SHA256SUMS.txt`：28/28 项通过。清单文件自身 SHA-256 为 `56522eaa220ac9a3c7ef87c9947f3a6a68ad4e35ffe699a7c4e399ada87866e1`。
- 使用的 `preflight_owner.py` SHA-256 为 `03234c8ab4dbb12d933d6c153d89e3be3af735d2650f00e11ab64516d4f678f3`。原始工具 JSON 见同目录的 [`D_PREFLIGHT_OWNER_RAW.json`](D_PREFLIGHT_OWNER_RAW.json)，未改写工具输出。

## 挂载、工作根与持久性

`findmnt -T /workspace` 返回 `TARGET=/workspace`、`SOURCE=overlay[/workspace]`、`FSTYPE=overlay`。挂载信息显示 Kata 容器 overlay；未发现单独挂载的持久工作卷。`/tmp` 是 tmpfs，不适合保存轨迹或 checkpoint。可见块设备没有独立挂载到工作目录；provider allocation、volume quota、跨正常会话重启持久性及恢复策略均无法从本机证明。`quota` 工具不可用，故 provider/用户级配额为 `UNVERIFIED`。没有访问 A/B 或其他云端目录。

当前候选路径均不存在，且没有将 overlay 冒充为通过验收的工作根：

- `sampler_work_root` 提案：`/workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot`；未创建、未设置。
- `DFT_WORK_ROOT`：沿用 C 固定工作根提案 `/workspace/dft-work/cloud-d/ag-ti-si-c-b1`；未创建、未设置。

这两个路径目前都落在同一未证实持久性的 overlay 上。应先提供已确认的持久卷，再将本机专属子目录设在该卷上；在此之前不创建目录、不把路径传给采样或 DFT runner。

## 资源与运行环境

在 `/workspace` 测量到可用 `28,797,059,072` 字节、总量 `33,770,192,896` 字节、可用 inode `2,038,540`。这是 overlay 的可见容量，不是经 provider 确认的卷配额。

CPU 可见数/affinity 为 5；`cpu.max=400000 100000`，有效 CPU quota 为 4 核。cgroup `memory.max=34,359,738,368` 字节、`memory.current=4,137,349,120` 字节；按可见限制计算，可用内存为 `30,222,389,248` 字节。系统 `MemAvailable` 为 `34,189,381,632` 字节；资源报告采用受 cgroup 限制后的较小值。更外层 provider 限制不可见。

未发现正在运行的 GPAW、LAMMPS、MACE、VASP 或 Quantum ESPRESSO 计算进程。此环境未提供调度队列读取接口；没有启动计算任务。

`/workspace/onboarding/runtime.sh` SHA-256 为 `d8b6d0c05020eb9baf270ea17ef67c0f8c043d5b6d5f5eac692eba007969e926`，`bash -n` 通过。未 source 或修改该脚本。现有配置要求 `GPAW_MPI_BACKEND=cgpaw`；`gpaw_data.datapath()` 解析为 `/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups`，与 `GPAW_SETUP_PATH` 一致。OpenMPI 为 5.0.7；设置该 backend 后，四 rank `gpaw.mpi.world.size == 4` 探针通过。

| 组件 | 实测版本/状态 |
|---|---|
| Debian | 13（现有环境） |
| Python | 3.12.14 |
| GPAW | 26.7.0 |
| ASE | 3.29.0 |
| gpaw-data | 1.2.1 |
| NumPy / SciPy | 2.5.3 / 1.18.1 |
| OpenMPI | 5.0.7 |

四元素 PAW 文件与固定项目 SHA 匹配：Ag `bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97`；Ti `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f`；Si `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0`；C `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983`。

## 容量预算与 checkpoint 保留

**Ag64 采样：**C 的提案对每条轨迹设 512 MiB 档案上限（含 restart bundle），D 有两项任务，合计上限为 1,073,741,824 字节；要求开跑前至少保留 2,147,483,648 字节。当前 overlay 可用空间在字节层面超过该预算。但因为没有已验证的持久根，`SAMPLING_STORAGE_READY=NO`。此处是容量比较，不是运行放行。

**GPAW checkpoint：**C 的资源说明要求预留至少 5,000,000,000 字节风险空间，并核算 current、previous、原子临时副本以及已保留 mesh 文件。当前 overlay 可用量高于这个单独的名义 reserve，但实际 64 Ag 冻结几何和 Gamma/shifted 2×2×2 checkpoint 大小尚未测量，持久根也未建立，因此 `DFT_STORAGE_READY=NO`。

保留台账状态：`NOT_STARTED_NO_DFT_ROOT`。两个提案根都不存在；D 本轮尚未生成 pilot checkpoint，也未删除、搬移或压缩任何数据。没有可声称通过的 DFT 保留台账。正式运行前须按 C 的 current+previous+atomic temporary 规则登记实测字节及哈希。

固定提交中 `C_RELEASE_CONFIG.json` 的 `storage_release=false`、`warmup_release=false`、D 的 `C_host_release=false` 保持原值。本报告不修改放行字段。未安装软件、下载势文件、生成轨迹或启动采样/DFT/训练。
