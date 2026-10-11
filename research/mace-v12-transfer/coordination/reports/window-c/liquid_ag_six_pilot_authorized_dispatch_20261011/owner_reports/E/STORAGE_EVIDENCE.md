# E 持久工作区与资源预检

采集时间：2026-10-11 09:42:38 Asia/Shanghai（2026-10-11T01:42:38Z）。

## 身份与证据来源

- `OWNER=E`：本环境存在 `/workspace/.setup/e/runtime.sh`；上传的 `START_E.md` 与固定提交中的版本逐字节相同，SHA-256 为 `6fd5e04beb69d4bdbfff21958df2d5f88b95d44ebb5f43a951d4b1e2b644f0ba`。这属于 E 自报身份依据，未获得 C 的远程主机认证。
- 固定提交：`2959fd3e2396a7319bd1f5c7f4fec27ce442abe0`。执行包 `SHA256SUMS.txt` 中 28 项全部通过；本报告目录另有独立的 `SHA256SUMS.txt`。
- `preflight_owner.py` SHA-256：`03234c8ab4dbb12d933d6c153d89e3be3af735d2650f00e11ab64516d4f678f3`。
- 原始 helper 输出见同目录 `RESOURCE_PREFLIGHT.json`。`identity_grade=OWNER_DECLARED; NOT_C_REMOTE_VERIFIED`，且 helper 明确记录持久性、资源隔离、provider quota 均未验证。

## 运行时与软件

- `/workspace/.setup/e/runtime.sh` SHA-256：`e32d19eef429939bcb3d0aa2d365d063b6dca16b4b4052b02da226fd7247e43b`；`bash -n` 通过，source 后解析出的 `GPAW_SETUP_PATH=/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups`，与 `gpaw_data.datapath()` 一致；`GPAW_MPI_BACKEND=cgpaw`。
- Python 3.12.14、GPAW 26.7.0、ASE 3.29.0、gpaw-data 1.2.1、NumPy 2.5.3、SciPy 1.18.1；`pip check` 无依赖错误。
- OpenMPI launcher 为 5.0.7。`gpaw info` 报告 MPI、ScaLAPACK、LibXC 5.2.3 可用；`_gpaw` 动态链接到当前 GPAW 环境中的 `libmpi.so.40`、`libscalapack-openmpi.so.2.2`、OpenBLAS、BLAS/LAPACK 与 LibXC，`ldd` 未见缺失库。只运行了四 rank GPAW MPI communicator 启动探针，四个 rank 均报告 `size=4`；未构造原子或计算能量/力。
- Ag/Ti/Si/C PBE setup SHA-256 均与项目基准一致：
  - Ag `bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97`
  - Ti `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f`
  - Si `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0`
  - C `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983`

## CPU、内存、磁盘与任务

| 项目 | 实测值 |
| --- | ---: |
| CPU affinity / cpuset | 5 个可见 CPU（0–4） |
| cgroup CPU quota | `400000 / 100000`，有效 4 CPU |
| cgroup memory limit | 34,359,738,368 B |
| cgroup memory.current | 3,848,351,744 B |
| 可用内存（取宿主 MemAvailable 与 cgroup 剩余量的较小值） | 30,511,386,624 B |
| `/workspace` 可用空间 | 28,798,726,144 B |
| `/workspace` 可用 inode | 2,039,455 |
| `/workspace` 文件系统容量 | 33,770,192,896 B |

`/workspace` 是可写 overlay 文件系统。`/workspace/scratch` 和 `/workspace/shared` 与它同一文件系统；`/tmp` 是 tmpfs；`/mnt/e-work`、`/data`、`/scratch` 均不存在或未挂载。provider 持久性策略、卷配额及卷级 inode 配额均没有可核验资料；overlay 视图不能证明正常会话重启后保留，也不能证明与 A/B、D/F 的底层资源隔离。没有执行重启测试。当前 PID namespace 未观察到 GPAW、LAMMPS、MD 或 DFT 作业进程。

因此，本次不在 `/workspace`、`scratch`、`shared` 或 `/tmp` 创建工作目录。候选独立卷需由 provider 挂载为 E 专属 `/mnt/e-work` 并提供持久性、配额和隔离依据。目录建议为：

- `sampler_work_root=/mnt/e-work/sampling`
- `DFT_WORK_ROOT=/mnt/e-work/dft`

两条路径目前均不存在，未标 PASS、未创建。现有 runtime 的 `e_configure_work_roots()` 会要求 `/mnt/e-work` 为单独可写挂载、相对 `/workspace` 为不同文件系统且至少有 30,000,000,000 B 空闲；当前挂载缺失，函数不会设置工作根。该函数当前把 DFT 根设为卷根，实际启用上述 `/dft` 子目录前还需由维护者对齐变量。

## 容量预算（采样与 GPAW checkpoint 分开）

**Ag64 采样：** 固定提交中的 E 配置仍标 `proposal_only_parameters=true`、`storage_release=false`、`warmup_release=false`。其中每条轨迹归档上限为 536,870,912 B；E 的两条轨迹合计上限 1,073,741,824 B，另有两份 33,554,432 B checkpoint reserve，合计约 1,140,850,688 B；配置的最低空闲线为 2,147,483,648 B。当前 overlay 的字节数在算术上高于这些提议门槛，但没有持久性、独立性或卷配额证明，所以 `SAMPLING_STORAGE_READY=NO`。没有启动采样。

**GPAW checkpoint：** 当前余量比既有 30,000,000,000 B 工作盘建议少 1,201,273,856 B。dispatch 的 storage note 要求至少保留 5,000,000,000 B 风险余量；因此本 overlay 上名义可用于新增数据的上限约为 23,798,726,144 B，但并非专属配额。实际 Ag64 的 bands、IBZ/PW/FFT 维度、checkpoint 字节数、mesh 数、current/previous/temporary 副本及累计保留量都未实测，因为坐标和 DFT 尚未放行。C 的 [磁盘资源预算说明](../../../five_cloud_preflight_and_liquid_pilot_readiness_20261011/DISK_AND_RESOURCE_BUDGET.md) 中 Gamma/2×2×2 保留情形约 16.773–23.402 GB、较密 3×3×3 情形约 31.032–50.422 GB，是含未验证假设的估算，不是 E 的 Ag64 实测需求。故 `DFT_STORAGE_READY=NO`，需在结构冻结后用真实 checkpoint 预算和专属持久卷重新验收。

## 状态与边界

```text
OWNER=E
RUNTIME_PASS=YES
PERSISTENCE_VERIFIED=NO
SAMPLING_STORAGE_READY=NO
DFT_STORAGE_READY=NO
RESOURCE_READY=NO
```

本轮未安装/升级软件，未创建采样/DFT 工作目录、结构或轨迹，未启动采样、DFT、训练或 TTM。未改动 `C_RELEASE_CONFIG.json`；其 SHA-256 仍为 `57dc948d45daa3c51f7b3e5491a19e4f407546806531c3afae14a221b054ffc0`，`storage_release`、`warmup_release` 和 `continuation_release` 均保持 `false`。E 保持 HOLD。解锁资源验收所需的外部事项是提供 E 独享的持久挂载、配额/恢复策略以及与 A/B、D/F 隔离的 provider 证据；完成后再用真实新挂载重跑预检。
