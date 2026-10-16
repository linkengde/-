# H：48 原子 Ti₃SiC₂ 静态小位移 DFT 单点预检修订

日期：2026-10-16。修订对象是 C 审核提交 `25bb74065797d516f1e8c8bbf9b68a4b558aefef` 中的意见，原始 H 预检报告保留不改。

```text
C 对原报告的结论      = NEEDS_REVISION
本修订状态             = 已提交 C 复审
C 对本修订的复审       = PENDING
DATA_ROLE              = HOLD
ENV_READY              = YES
RESOURCE_READY         = NO
DFT_EXECUTION_READY    = NO
execution_authorized   = false
```

本轮只做离线输入、k 点、软件/PAW/MPI 与资源预检。没有启动 DFT、MD、AIMD、TTM、松弛或训练，没有修改结构、晶胞、原子顺序、历史计算或 checkpoint，也没有创建 DFT 输出目录。

## 1. 固定输入与数据角色

- 唯一候选：`TSC48_static_perturb_std0p025_seed2026101603_UNLABELED.extxyz`
- 输入来源提交：`cc65a3e6b75be28c06fda8608b4064662d714117`
- 输入 SHA-256：`da3820c477023e0bb423f87550ed82ae0c2a571f53d294c055f78d7d052d5e8b`
- 48 原子；母体输入 SHA-256：`157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4`
- 候选与母体晶胞、PBC、组成及原子顺序未改。候选仍是 `TSC_BULK_PARENT_OR_CONTROL_FAMILY` 的静态扰动，保持 `DATA_ROLE=HOLD`；它不是高温结构，也不得独立划为 DEV/TEST。

原预检报告路径为 `../H_TSC48_STATIC_PERTURB_SINGLEPOINT_PREFLIGHT.md`，SHA-256 `eec7bb2b0af188a9cd34f8199c234ccca361a53af1e533a23d7f06c501d31aa8`（与本分支 HEAD 中的版本一致）。本修订另存于新目录，没有覆盖它。候选文件及既有计算结果均未改。

## 2. k 点定义已纠正

| 项目 | 原 H 报告 | 本修订候选 |
|---|---|---|
| 网格大小 | 10×10×4 | 保留 10×10×4 |
| 网格定义 | 错写为 Γ 居中 | 标准偶数 Monkhorst–Pack，与母体实际实现一致 |
| GPAW 参数约定 | 将 shift=0 误当作 Γ 居中 | 精确沿用母体 runner 的 `kpts=tuple(record['method']['kpts'])`；不传 `gamma=True`，不传额外 offset |
| 实际 shift | 报告称 Γ 居中 | `[0,0,0]`；偶数 MP 网格不包含 Γ 点 |
| 适用结论 | 曾用母体收敛结果支持错误定义 | 母体网格收敛只支持该母体的网格比较；候选扰动几何自身的 k 点收敛仍为 `UNVERIFIED` |

**母体实现核对：**母体输入为 48 原子 `Ti3SiC2_baseline_PROPOSAL.extxyz`，SHA-256 `157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4`。母体 runner SHA-256 为 `b741d3b10193ad3bc7d63ae4bc20364f373f5e0bbe75bc9de32edfcdb0f858d8`，实际将整数三元组传给 GPAW。母体原始日志 SHA-256 为 `bb55d1933ad649e5a5c2cb1182108c8b5e3a8fb42cd0eeb3d95b91fd690a3dce`，记录 `MonkhorstPack(size=[10,10,4])`、MP shift `[0,0,0]`、400 个 BZ 点和 200 个 IBZ 点。日志中可见的首 10 行及最后一行 k 点坐标均与离线生成的母体 IBZ 一致，例如首行为 `(0.05,-0.45,-0.375)`、末行为 `(0.45,0.45,0.375)`；日志对中间行以省略号显示。

**离线实现核对：**读取 GPAW 26.7.0 `gpaw.dft.MonkhorstPack.build()` 和 GPAW/ASE 的 `kpts2sizeandoffsets` 实现。该接口默认 `gamma=None`；只有显式提供 `gamma` 才计算非零 offset。用母体 runner 的整数网格约定离线生成完整坐标，使用的标准 MP 坐标为 `(index + 0.5) / size - 0.5`。母体与候选的 400 个 BZ 坐标数组逐项完全相同；两者离线约化后也各有 200 个 IBZ 点。母体有 2 个空间对称操作，扰动候选有 1 个；候选的 200 点计数利用时间反演对称约化。这些只是坐标/对称性与资源预检，不是 DFT 计算或候选收敛证据。

完整 400 点候选网格表：[`candidate_bz_10x10x4_full.csv`](candidate_bz_10x10x4_full.csv)，SHA-256 `a59b5caefcd50109b408f8b4112e452f8abb2966c96bbbf69991183f3a6c16e0`。母体日志可见坐标摘录：[`parent_log_visible_ibz_excerpt.csv`](parent_log_visible_ibz_excerpt.csv)，SHA-256 `fa6f8d1a371bfa68f2dfe133504c65e8ab730794389de98b385257ac5c00f5ab`。机器可读审计：[`kpoint_grid_audit.json`](kpoint_grid_audit.json)。

## 3. 单点配置锁定与哈希

候选输入仍固定为原 SHA，不松弛，不改变晶胞或坐标。供 C 审核的唯一候选配置如下；这是静态配置记录，不是启动脚本：

| 参数 | 锁定值 |
|---|---|
| XC / 基组 | PBE / 平面波 |
| 截断能 | 500 eV |
| 占据 | Fermi–Dirac，0.10 eV |
| 自旋 | 非自旋极化 |
| k 点 | 标准 MP `(10,10,4)`；gamma 未指定；offset `[0,0,0]` |
| 本征态求解器 | PPCG，显式 `niter=2` |
| Mixer | `Mixer(0.05, 8, 100)`，即 Pulay backend，beta=0.05、nmaxold=8、weight=100 |
| SCF 收敛阈值 | energy、density、eigenstates 均 `1e-5` |
| 最大 SCF 迭代 | 160 |
| MPI | 4 ranks；每 rank 的 OMP/BLAS 线程数建议为 1 |

锁定配置文件：[`calculation_config_preflight.json`](calculation_config_preflight.json)，SHA-256 `ca6f365e8b0b8ec847e329d0312ebb51381b21e996050949494f2d3e0669dfe7`。

当前环境版本与构建文件核验：Python 3.12.14、GPAW 26.7.0、ASE 3.29.0、gpaw-data 1.2.1、OpenMPI 5.0.7。`GPAW_SETUP_PATH=/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups`，与 `gpaw_data.datapath()` 相同。运行脚本 `/workspace/.setup/runtime.sh` SHA-256 为 `e5e36ce44bb9e2e0f0561aadde291b393e6f2ac49605f2b82ac5f86baf0c74b1`。

- GPAW `__init__.py` SHA-256：`fb370df2212354f5af291761935988c9218602b7f3418d5f89e6a1cb66995005`
- ASE `__init__.py` SHA-256：`92ce3cc614521afe055e962ccb4254e8fbf342f3ef86111fa4abc0b017a9963e`
- GPAW `_gpaw` 扩展 SHA-256：`417bab86c51a5ccde9e48f8f863d76efaea5f3bcffc2dd175ce83dd66af1785a`
- `mpiexec` 5.0.7 可执行文件 SHA-256：`96ab74bcbe60c8b6a750cd3b726ce6b2360991bcbc2a16635525a1979fc252ff`
- `_gpaw` 当前链接到 `libmpi.so.40`，实际库 SHA-256：`dcca32f78fe6bda5ab10ca9019a3129e1869ffe6bb7a6f5461207563d2cea446`；OpenMPI ScaLAPACK `libscalapack-openmpi.so.2.2` SHA-256：`dc09400498e5806e6ef05e7278fab28820d6a74caa36ba5a804f7c96ff1e3939`
- Ti/Si/C PBE PAW SHA-256：Ti `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f`；Si `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0`；C `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983`。

所有版本、运行二进制/库哈希和 PAW 路径见 [`runtime_abi_hashes.json`](runtime_abi_hashes.json)。本轮用 `mpiexec -n 4` 运行 GPAW MPI 环境探针；4 个 rank 的 Python/GPAW/ASE/gpaw-data 版本、runtime 路径和 Ti/Si/C PAW 哈希一致，公共状态 SHA-256 为 `07c29cc90903c6c6ffdb1d64b467e8e8b2926dc91ac486b82121a81e2ff9adb9`。探针未创建 GPAW calculator，也未启动 DFT。MPI 环境一致性：`PASS`；它不等同于本候选的 4-rank DFT 运行测试。

## 4. 磁盘、内存、checkpoint 与墙钟预算

资源快照记录在 [`resource_budget.json`](resource_budget.json)，测量时间 `2026-10-16T10:08:25Z` 附近：

| 项目 | 只读证据 | 判断 |
|---|---:|---|
| CPU 配额 | `cpu.max=400000 100000`，等效 4 核；cpuset `0-4` | 4 ranks 与配额相符 |
| 内存 | 上限 34,359,738,368 bytes（32 GiB）；采样时已用 4,168,507,392 bytes，限额内余量 30,191,230,976 bytes | 母体同规模日志可供参考；本扰动任务峰值未实测 |
| `/workspace` 可用空间 | 28,827,836,416 bytes（28.83 GB / 26.85 GiB） | 比 30 GB 建议余量少 1,172,163,584 bytes；30 GB 是安全建议，不是 GPAW 固定最低要求 |
| scratch | `/workspace/scratch` 目前为空，与 `/workspace` 是同一 overlay 文件系统 | 单点 scratch 峰值未知；不能算作额外磁盘 |
| checkpoint 持久性 | `/workspace` 是 overlay；`/tmp` 是 tmpfs；没有核实到独立持久卷或跨环境重启保留证据 | 持久性、恢复可用性 `UNVERIFIED` |
| 外部墙钟 | shell `ulimit -t` 为 `unlimited`；当前环境没有可读的平台级作业超时证据 | 不代表会话墙钟；至少 5 小时配额 `UNVERIFIED` |

**磁盘预算估算（不是候选实测）：**历史 48 原子 8×8×4 的有效 `state.gpw` 为 7,361,378,540 bytes，128 个 IBZ 点，历史文件 SHA-256 `3a563ea6def61886af7d6be7547b591fad510e6778daff5b75c678f61d0563dc`。按当前候选离线得到的 200 个 IBZ 点线性放大，单份 `.gpw` 估算约 11,502,153,969 bytes。为保留旧的最后一个有效 checkpoint 并暂存新 checkpoint，写入峰值约 23,004,307,938 bytes；扣除后约余 5,823,528,478 bytes 供 scratch、日志、临时文件及其他占用。母体 10×10×4 的归档结果文件合计为 31,681 bytes，但这不包含大型 checkpoint。母体日志报告 wavefunction 10,907,008,000 bytes、projectors 3,388,134,400 bytes。扰动输入可能改变执行行为；候选 checkpoint 与 scratch 峰值没有实测或严格上界，因此不能据此证明预算一定足够，也暂不申请低于 30 GB 建议值的例外。

**checkpoint/恢复方案建议（未配置、未测试）：**先由平台提供并确认 H 专属持久工作卷，再把不可变输入与上述配置哈希写入该卷的独立作业目录。按母体 runner 的 20 步节奏候选保存到同卷临时名 `state.next.gpw`；只有所有 rank 写入完成、文件检查通过并完成落盘后，才原子替换 `state.gpw`，确保写入失败时旧的最后有效 checkpoint 仍在。恢复时先复核输入、配置、PAW 哈希，再通过 GPAW `restart()` 从 `state.gpw` 读入并继续 SCF；该恢复流程尚未在此环境验证。发生 ENOSPC、挂载丢失或 checkpoint 不完整时，停止任务并标记 `FAILED/INCOMPLETE`，保留现有文件供审核，不从临时/损坏文件恢复，不发布能量或力标签，不删除历史科研数据。大型 `.gpw` 不应推送到 GitHub；GitHub 仅用于小型报告与校验清单。

因此当前仍为 `RESOURCE_READY=NO`。30 GB 例外须等 checkpoint 写入峰值、scratch 上界、持久卷与恢复方式、至少 5 小时的平台墙钟证据均明确后，再由 C 审核。

## 5. 尚未解决事项与需要 C 再审的项目

1. C 需确认修订后的标准 MP 定义、离线 k 点列表/哈希及锁定配置；本修订提交前，原 `NEEDS_REVISION` 结论仍有效。
2. 候选结构自己的 k 点收敛仍未验证；母体网格收敛不能替代此验证。本次只提出首个候选网格，不声称误差已收敛。
3. 需要平台级外部墙钟上限证据（建议至少 5 小时），以及 H 工作区/卷跨重启持久性的证据。
4. 需要限制并测量 scratch 和 checkpoint 临时写入峰值，确认当前剩余空间足够；必要时由 C 审核 30 GB 安全建议的例外。
5. 需 C 确认独立工作根和 checkpoint/恢复/存储故障方案；本报告未创建工作根或 checkpoint。
6. 数据仍保持 `DATA_ROLE=HOLD`，归属 `TSC_BULK_PARENT_OR_CONTROL_FAMILY`；不得自行转入训练或作为独立 DEV/TEST。
7. 只有 C 完成科学与资源复审、资源正式放行，且用户对这一唯一单点另行明确授权后，才能准备执行。当前没有执行授权。

## 6. 停止状态

本修订只提交预检记录与坐标、配置、运行时和资源证据。`ENV_READY=YES`，但 `RESOURCE_READY=NO`、`DFT_EXECUTION_READY=NO`。任务在等待 C 复审；没有自动启动任何计算或下一项任务。
