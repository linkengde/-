# H：48 原子 Ti₃SiC₂ 静态小位移单点 DFT 预检

**状态：`READY_FOR_C_REVIEW_ONLY` · `RESOURCE_READY=NO` · `DFT_EXECUTION_READY=NO`**

本报告仅记录只读预检和一项候选单点方案。未启动 DFT、MD、AIMD、松弛、TTM 或训练；未修改输入结构或已有标签。用户要求的 C 审核和明确计算授权均未获得，任务保持 HOLD。

## 1. 输入结构核验

- 源提交：`cc65a3e6b75be28c06fda8608b4064662d714117`
- 源文件：`research/mace-v12-transfer/coordination/reports/chatgpt-nep-ace-dft-repair-20261016/NEW_DFT_GEOMETRY_CANDIDATES/structures/TSC48_static_perturb_std0p025_seed2026101603_UNLABELED.extxyz`
- 文件 SHA-256：`da3820c477023e0bb423f87550ed82ae0c2a571f53d294c055f78d7d052d5e8b`（与用户给定值一致）；大小 2,572 bytes。
- ASE 3.29.0 只读解析：48 原子，化学式 `C16Si8Ti24`，`parent_row_id` 连续为 1–48；记录中没有能量或力标签，`label_status=NO_DFT_ENERGY_OR_FORCES`、`role=DFT_INPUT_ONLY`。
- 晶胞（Å，行向量）：`[[6.115, 0, 0], [-3.0575, 5.2957453441, 0], [0, 0, 17.6235]]`；晶格长度 `[6.115, 6.115, 17.6235] Å`，夹角 `[90°, 90°, 120°]`；PBC=`[true,true,true]`。
- 使用周期性最小镜像计算的最近原子对为 Ti–C，距离 `2.0016108 Å`；ASE 0-based 索引为 5、23（输入行号 6、24）。
- 对照同提交中的母体 `Ti3SiC2_baseline_PROPOSAL.extxyz`：元素顺序和晶胞相同，母体 SHA-256 为 `157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4`；位移 RMS `0.0416083 Å`，最大单原子位移 `0.0722999 Å`。这些只是静态坐标差，不代表高温结构或热采样结果。

**结构检查：PASS。** SHA、原子数、组成、晶胞、PBC 和周期最近邻均通过基本输入核验；这不等于已经计算或验证该结构的能量、力或热力学稳定性。

## 2. Ti₃SiC₂ 体相 k 点证据与候选网格

仓库已有同一 48 原子 CIF 派生母体、同一 PBE/PW500 设置的固定几何网格比较。`8×8×4 → 10×10×4` 只改变 k 网格，两个端点的结构输入 SHA 相同且归档校验通过：

| 指标（10×10×4 减 8×8×4） | 母体实测差 | 项目基准 |
|---|---:|---:|
| Native 总能 | −0.03512 meV/atom | ≤ 2 meV/atom |
| Force-consistent 自由能 | −0.01669 meV/atom | ≤ 2 meV/atom |
| 全原子力向量 RMS | 0.004618 eV/Å | ≤ 0.01 eV/Å |
| 最大单原子力向量差 | 0.007105 eV/Å | — |

更粗的 `6×6×4 → 8×8×4` 力 RMS 为 `0.010850 eV/Å`，超过原始 0.01 eV/Å 门槛；`6×6×4 → 10×10×4` 为 `0.006258 eV/Å`，通过该母体上的差值门槛。该证据支持把 **Γ 居中 `10×10×4`** 作为本输入的一次候选单点网格；这是按同一 48 原子体相晶胞证据提出的设置，不是对扰动结构已收敛的声明。现有直接比较没有覆盖本扰动几何，因此其自身 k 点误差仍为 **UNVERIFIED**。此处没有采用 32 原子界面的 `1×10×8` 网格。

证据位置：

- [`k8x8x4_k10x10x4_comparison.json`](../../window-a/v19_Ti3SiC2_bulk_mesh_comparison/k8x8x4_k10x10x4_comparison.json)
- [`mesh_matrix_kz2_kz4.json`](../../window-a/v19_Ti3SiC2_bulk_mesh_comparison/mesh_matrix_kz2_kz4.json)
- 母体 10×10×4 归档：`research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/calculations/Ti3SiC2_baseline_k10x10x4_v19/`

母体 10×10×4 日志记录 `MonkhorstPack(size=[10,10,4])`、shift `[0,0,0]`（Γ 居中），SCF 34 步，4 个 MPI ranks，用时 `11,722.32 s`（约 3.26 h）。母体结果不能替代扰动结构的收敛验证。

## 3. 与已有 48 原子 DFT 几何的去重和父族检查

`DFT_MASTER_INVENTORY.csv` 中已有五个唯一、完整且 SCF 收敛的 48 原子 Ti₃SiC₂ 标注几何：

1. `Ti3SiC2_baseline_k10x10x2_v19`（同一母体几何有多个网格标签，不能按网格变体重复计作多个结构）；
2. `Ti3SiC2_Ti4f_id0003_z_p020A_k10x10x4_v19`；
3. `Ti3SiC2_Ti4f_id0003_z_m020A_k10x10x4_v19`；
4. `Ti3SiC2_C4f_id0009_z_p020A_k10x10x4_v19`；
5. `Ti3SiC2_C4f_id0009_z_m020A_k10x10x4_v19`。

五者的 `parent_configuration_family_id` 均为 `TSC_BULK_PARENT_OR_CONTROL_FAMILY`，母体为 `V19_Ti3SiC2_48ATOM_BULK_CONTROL_FAMILY`，标签状态为 `VALID_DFT_LABEL`，能量/力完整且 SCF 收敛。按元素顺序、晶胞和坐标比较，本候选与五个已有几何均非 `1e-6 Å` 精度下的精确重复；它是相同母体族的另一个静态扰动点。

因此，若 C 后续同意纳入训练，应沿用该 Ti₃SiC₂ 母体族标识，并在数据拆分时把全族放在同一分区；不得把候选点当成独立 DEV/TEST 样本，也不得把母体不同 k 网格标签重复计数。清单中的训练/验证/测试角色尚未冻结，需由 C 决定最终角色。

清单位置：[`DFT_MASTER_INVENTORY.csv`](../../window-e/nep_ace_dft_qualification_audit_20261015/DFT_MASTER_INVENTORY.csv)。

## 4. PBE、PAW、MPI 和运行环境

当前只读环境检查及仓库归档证据：

| 项目 | 核验事实 |
|---|---|
| Python / GPAW / ASE | Python 3.12.14 / GPAW 26.7.0 / ASE 3.29.0 |
| GPAW data / 数值库 | gpaw-data 1.2.1；NumPy 2.5.3；SciPy 1.18.1；LibXC 5.2.3；ScaLAPACK 可用 |
| 方法基准 | PBE；平面波截断 500 eV；Fermi-Dirac 展宽 0.10 eV；非自旋极化 |
| SCF 基准 | 能量、密度、特征态阈值各 `1e-5`；Pulay `beta=0.05, nmaxold=8, weight=100`；`maxiter=160` |
| PAW 包 | Ag/Ti/Si/C PBE 文件 SHA-256 与项目基准一致：Ag `bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97`；Ti `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f`；Si `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0`；C `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983` |
| MPI / ABI | OpenMPI 5.0.7；`_gpaw` 链接到 `libmpi.so.40` 和 `libscalapack-openmpi.so.2.2`，当前运行库均已解析。H 既有 4-rank MPI 测试及 GPAW 任务通过，rank 0–3 均报告 size=4；本次未重跑测试。 |
| runtime.sh | `/workspace/.setup/runtime.sh` SHA-256 `e5e36ce44bb9e2e0f0561aadde291b393e6f2ac49605f2b82ac5f86baf0c74b1`；source 后 `GPAW_SETUP_PATH=/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups`，与 `gpaw_data.datapath()` 相同。 |

环境可用性结论：`ENV_READY=YES`，通过已验证的绝对 Python 路径和 runtime 脚本使用。脚本没有把 venv 的 `bin` 加入默认 `PATH`，因此裸 `gpaw` 命令不可用；候选启动方案应使用 `/workspace/.venvs/gpaw-mpi/bin/python` 或显式 venv `bin/gpaw`，不改全局配置。已归档的 4-rank GPAW 任务证明 MPI 可用，但不代表本次已执行任何计算。

环境证据：`research/mace-v12-transfer/coordination/reports/chatgpt-nep-ace-dft-repair-20261016/CLOUD_DFT_RESULTS_20261016/H_Ti_gap_32_k1x10x8/preflight/{VERSIONS_AND_RUNTIME.txt,MPI_4RANK_TEST.txt,GPAW_MPI_LINKAGE.txt}`。

## 5. 资源和时间限制

| 资源 | 只读观察 | 评估 |
|---|---:|---|
| CPU | cgroup 配额 4 CPU；cpuset 显示 `0-4`，`cpu.max=400000/100000` | 4-rank 作业可用；CPU 亲和集合比配额多一个逻辑 CPU，不等于 5 核吞吐 |
| 内存 | 上限 `34,359,738,368` bytes（32 GiB）；检查时已用 `4,122,353,664` bytes，限额余量约 28.16 GiB | 尚无本扰动任务峰值 RSS 记录，峰值需求为 UNVERIFIED |
| `/workspace` 可用磁盘 | `28,828,061,696` bytes（28.83 GB / 26.84 GiB） | 低于项目建议的 30 GB，短缺约 1.17 GB；不删除历史文件扩容 |
| 持久化 / checkpoint 盘 | 未发现独立持久盘挂载；`/workspace` 为 overlay，`/tmp` 为 tmpfs（约 17.30 GB 可用） | `/tmp` 不作为持久 checkpoint 位置；`/workspace` 的跨任务持久性仍需平台侧确认 |
| 进程 CPU 时间 | `ulimit -t` 为 unlimited | 仅说明该 shell 的 CPU 时间限制，不说明云端作业/会话硬超时 |
| 外部墙钟时限 | 无可验证的调度器/任务时限记录 | **UNVERIFIED** |

可比母体的 10×10×4 实测约 3.26 h；8×8×4 曾用约 4.23 h。它们不是对本扰动任务的时长保证。建议 C 只在确认单作业墙钟配额不少于 5 h 后考虑放行，并保留余量。本任务 checkpoint 峰值空间及 `/workspace` 的跨任务持久性均未确认。

因此 `RESOURCE_READY=NO`：磁盘未达到 30 GB 项目建议线，外部墙钟时限和 checkpoint 持久化属性也未核实。

## 6. 唯一候选单点方案（尚未授权、未执行）

| 项目 | 候选设置 |
|---|---|
| 唯一输入 | `TSC48_static_perturb_std0p025_seed2026101603_UNLABELED.extxyz`，固定使用上述 SHA；不重写、不松弛、不改变晶胞或坐标 |
| 作业标识 | `H_TSC48_static_perturb_std0p025_seed2026101603_gamma10x10x4` |
| 网格 | Γ 居中 `10×10×4`；依据同晶胞母体收敛证据提出，候选本身网格误差仍待评估 |
| 电子结构 | PBE / PW 500 eV / Fermi-Dirac 0.10 eV / 非自旋极化 |
| SCF | 能量、密度、特征态阈值均 `1e-5`；Pulay `0.05/8/100`；最多 160 步 |
| 并行 | 4 个 OpenMPI ranks；沿用已验证的 H runtime 和 PAW 文件 |
| 输出隔离 | 仅在 H 专属、已确认持久且空间合格的 `DFT_WORK_ROOT/TSC48_static_perturb_std0p025_seed2026101603/gamma_10x10x4/` 新目录写入；不得覆盖现有数据。该根目录尚未创建，持久挂载尚待核实。 |

**放行前必须通过：** C 确认数据角色/父族归属并审核候选网格；用户给出该唯一单点的明确计算授权；可用磁盘达到项目门槛且 checkpoint 挂载持久；确认至少 5 h 外部墙钟配额；再次核对输入 SHA、PAW SHA、GPAW_SETUP_PATH、4-rank 启动和空输出目录。未满足前保持 HOLD。

`DFT_EXECUTION_READY=NO`。本报告不包含也不调用计算启动命令。

## 7. C 审核状态复查

用户要求等待 C 审核。只读检查当前可访问的 C 相关远端分支时，分支头分别为 `cloud-c-tisc-dft-extension-plan-20261014`（2026-10-14）、`cloud-c-ag-tisc-cross-dft-review-20261015`（2026-10-15）及 `cloud-c-five-cloud-dft-plan`（2026-10-10）；这些分支和 `main` 均未检索到本候选文件名、作业标识或对应预检批准记录。因此 C 对本候选的审查状态仍是 **PENDING / 未找到记录**。本次“继续”未被解释为越过该审核门槛的单点计算授权。

## 8. 总结状态

| 检查项 | 状态 | 说明 |
|---|---|---|
| 输入 SHA / 结构基本有效性 | PASS | 48 原子、PBC、晶胞及最近邻已核验 |
| 已有 Ti₃SiC₂ 母体网格依据 | PASS（限母体） | 10×10×4 在原始能量/力差门槛内；不证明候选扰动点已收敛 |
| 现有几何去重 / 父族归属 | PASS（候选不重复；同一父族） | 需 C 冻结训练/评估角色 |
| 软件 / PAW / MPI / runtime | PASS（绝对路径） | runtime 未配置 venv `bin` 默认 PATH |
| 资源 / 墙钟 / checkpoint | NOT READY | 磁盘低于 30 GB；墙钟及跨任务持久挂载 UNVERIFIED |
| DFT 执行 | **HOLD / NO** | 等待 C 审核及用户明确授权 |

静态小位移候选仅是一份固定几何输入；本检查没有、也不能将其描述为高温 Ti₃SiC₂ 结构。
