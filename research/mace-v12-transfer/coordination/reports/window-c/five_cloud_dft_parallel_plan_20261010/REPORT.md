# 五云端 DFT 并行计划：Ag–Ti₃SiC₂ 统一势函数

日期：2026-10-10
状态：**仅规划与只读预检设计；未放行计算。** 本交付不生成结构、不启动 DFT/MD/AIMD/TTM/MACE 训练，不改势函数、数据角色、封存数据或 A/B 任务。

## 结论

把原 30 个拟议单点按 D/E/F/G/H 五个独立云端登记为 18 个液态 Ag 候选、6 个 Ti₃SiC₂ 骨架候选、6 个固态 Ag 候选。每个候选有唯一的提案 task ID 和 parent-family ID；这些是规划编号，不是已启动作业号。完整 30 行见 `five_cloud_task_matrix.csv`。

现阶段**没有条件启动六个 1250 K 液态 Ag pilot**：坐标和轨迹均未生成，1250 K 液态 Ag 参考密度的精确文献值/方程尚未核验，候选采样器和许可也未冻结；同时当前没有 D/E/F 主机资源或运行队列的可见连接。六点的矩阵位置已由 D/E/F 各承担两个，G/H 不参与首轮 pilot。必须先完成来源、密度、主机报告、实际结构及液态判据审核，并由用户另行批准结构生成及具体 DFT pilot。

D–H 的独立资源关系来自用户确认；C 当前没有访问这五台主机的 SSH/云实例工具。当前 `origin/main` 快照只有 A/B 协作记录，没有 D–H 任务卡；这不能证明 D/E/F 没有本地运行任务。开始重新分配前，各云端 owner 必须只读报告正在运行、已排队和未启动的工作。任何已运行任务和 checkpoint 保持原样，只安排尚未启动的构型。

云端 C 的分析环境已按仓库锁定版本在隔离的 `/tmp` 虚拟环境中准备：Python 3.12.14、PyTorch 2.5.1+cpu、mace-torch 0.3.16、ASE 3.29.0、NumPy 1.26.4、SciPy 1.17.1；`pip check` 与 `mace_run_train --help` 均通过。没有运行训练或模型推理，也没有安装 GPAW。此环境结果只说明 C 协调执行器可用，不说明 D–H 已安装或兼容。

## 任务分配

| 云端 | 液态 Ag（18） | Ti₃SiC₂（6） | 固态 Ag（6） | pilot 阶段 | 拟议总量 | 初估串行 wall time* |
|---|---:|---:|---:|---|---:|---:|
| D | 4 | 1 | 1 | 1250 K 两点 | 6 | 17–42 h |
| E | 4 | 1 | 1 | 1250 K 两点 | 6 | 17–42 h |
| F | 4 | 1 | 1 | 1250 K 两点 | 6 | 17–42 h |
| G | 3 | 2 | 1 | 预留；先不跑液体 | 6 | 17–40 h |
| H | 3 | 1 | 2 | 预留；先不跑液体 | 6 | 16–38 h |
| **合计** | **18** | **6** | **6** | **D/E/F 六个液体候选** | **30** | **84–204 h** |


data-level 分配、parent、seed、温度、密度因子、候选来源与状态均在 CSV。六个 pilot 为 1250 K、0.97/1.00/1.03 × 两个 seed，每个 64 个 Ag 原子；配对分散到 D/E/F：

- D：0.97 × seed01；1.03 × seed02。
- E：0.97 × seed02；1.00 × seed01。
- F：1.00 × seed02；1.03 × seed01。

每个 seed 是独立采样 parent。其余 12 个液体候选（1050 K、1600 K 两温度 × 三密度因子 × 两 seed）暂列为 pilot 后候选；G/H 上的任务以及 D/E/F 的非 pilot 任务均保持未放行。温度专属液态密度来源需分别确认，不能默认 1250 K 密度适用于 1050/1600 K。

六个 Ti₃SiC₂ 候选拟来自两个新 parent families、每个三个骨架/模式结构；六个固态 Ag 候选拟来自两个新 parent families、每个三个小扰动/响应结构。精确超胞、原子数、坐标、来源 hash、seed 和体积均待结构方案冻结后登记。不能选用已经评分的父族充当独立验证，也不改变旧纯相数值控制标签的角色。

## 资源与预计时间

D–H 每台的 CPU、内存、剩余磁盘、OS、Python、GPAW、ASE、MPI ABI、BLAS/ScaLAPACK/LibXC、PAW 数据版本/许可/校验值和实际队列状态目前均为 **UNVERIFIED**，见 `resource_preflight.csv`。本计划不声称已经远程检查或安装这些软件。

负载估计沿用归档锚点：相近的 48–52 原子 GPAW 作业约 2.8–3.8 h；四 MPI ranks 下暂按液态 Ag 64 原子 3–8 h/例、Ti₃SiC₂ 48–96 原子 3–6 h/例、固态 Ag 32–64 原子 2–4 h/例。30 个最终单点合计约 84–204 个四-rank 串行 wall 小时，即 336–816 rank-hours。单云端估算为上表范围。按五台同时可用且每台一次只运行一个 4-rank 作业计算，理想 makespan 由最慢云端约 17–42 h 决定；这是算术估算，不是性能承诺，不包含安装、排队、候选生成、mesh pilot、异常重跑或硬件差异。

六点 pilot 的六个 DFT 单点约 18–48 h；同一 64 原子液体快照的 Γ/2×2×2 数值对照增加约 3–16 h，合计 21–64 个四-rank 串行 wall 小时（84–256 rank-hours）。按 D/E/F 同时各跑两点、D 另承担 1–2 个 mesh 对照，条件性并行 makespan 约 9–32 h。首例若超过 8 h，先停下评估成本与方法；每个获准单点建议上限 12 h 或 160 次 SCF iteration，触及任一上限即保存日志/checkpoint 并停批，不自动重启。

建议 pilot 前每台至少确认 4 个可用物理核、液相 pilot 至少 16 GiB 可用内存；若运行 80–160 原子后续作业，建议 32 GiB 以上；每个独立计算区建议至少 30 GB 可用空间。实际适用值须结合 GPAW/PW/k 网格、队列限制和预检输出确认。依据用户建议的 DFT 配置，4–8 核、32 GB 以上内存及 Linux 是较稳妥目标；GPU 对 GPAW CPU 路线非必需。

## 方法与标签统一

候选方法为 GPAW 26.7.0、ASE 3.29.0、PBE、PAW data 1.2.1、PW 500 eV、电子展宽 σ=0.10 eV、SCF 目标 1e-5。离子采样温度和电子展宽必须分开记录。每条结构需保存原始/最终结构哈希、元素/原子 ID 与顺序、cell/PBC、parent/seed/trajectory lineage、采样器及版本/许可、温度/密度、输入设置、k 网格、SCF 收敛、能量、力和结果哈希。能量单位 eV、力单位 eV/Å；用于监督的能量键与力须同一口径，热结构优先统一为 `free_energy` + `REF_forces`，另保留 native `energy`，不得把二者混用。

液体/共存盒先对同一个固定快照做 Γ 与 2×2×2 k 网格比较；若 ΔE > 2 meV/atom、force RMS > 0.01 eV/Å 或最大逐原子差 > 0.03 eV/Å，则尝试 3×3×3 或更密设置并重新审查。上述数值是拟议预算门，不是已验证精度。纯固相沿用同 parent 已收敛的倒格矢密度并复核超胞尺度，不能所有体系一律 Gamma。每个云端使用与 GPAW 编译 ABI 匹配的 MPI，不混合 OpenMPI/MPICH runtime；记录 GPAW 构建、BLAS/LAPACK、ScaLAPACK 和 LibXC 能力。PAW data 的具体版本、来源、许可、完整性校验以及可用于项目的条款必须逐台核实。

现有审计记录显示已归档的 15 个 Ag 和 10 个 Ti₃SiC₂ 数值控制标签仍保持 `numerical_pure_phase_control_not_training_or_independent_validation`。本轮不改变它们角色，也不复制进训练集。公开 DFT 界面记录需要按原报告的 hash/角色审核，不能因进入 30 点矩阵而自动变为训练标签。30 个新点均只是**未来训练候选**；全部父族须在输出标签可见前冻结。新独立开发/测试必须使用从未评分的新父族；这些 30 点不构成独立 test。封存测试保持不可见。

## 重复、父族与文件隔离

几何生成后先比较输入文件 SHA-256、元素顺序、cell/PBC、规范化周期结构指纹及父源/轨迹 lineage。液体任何同轨迹快照归在同一 parent；不同 seed 即不同采样 parent，但需确认不是同一 melt 轨迹续跑。纯相所有同一母晶体衍生扰动归成同一 parent。发现与旧归档或别云端任务同几何时，暂停该 task ID 并报协调，不静默复算/重标。

每台云端以自己的 branch、任务表、输入目录、checkpoint、输出、日志和 SHA 清单隔离保存。输入若需跨主机镜像，复制只读冻结输入包并核对 SHA；每台仍使用自己的 run/checkpoint 路径。不得把活跃 checkpoint 放进共享目录，或把大型 `.gpw` 文件提交 GitHub；只在通过核验后提交紧凑输出和 provenance。C 只维护协调清单，不改 A/B 工作区或任务卡。

## 阶段与停止条件

1. **主机与数据预检：** D–H owner 提交只读资源报告；D/E/F 额外标明已有/运行中/排队/未启动任务。任何现有进程、输入和 checkpoint 原样保留。没有真实坐标、来源、density basis、parent/seed 或方法/许可信息时，不创建 DFT 输入。
2. **候选液体准备：** 1250 K 熔点参考约 1234.93 K；精确 `rho_liquid(1250 K)` 尚未在可读公开来源中核实。液态坐标尚不存在。选定经来源/许可核验的 sampler 后，须经另行批准生成候选轨迹；MACE 或 Ag-only EAM 仅可供结构提议，不能把势能/力作为 DFT 标签。使用温度/体积/能量/应力平稳、RDF 与结构因子无持久 fcc Bragg、CN 分布、Q6/CNA、非平台 MSD 和 seed 间 RDF/CN/SOAP 差异证明是真液体。64 原子周期盒只为局部标签，不验证熔点/凝固。
3. **pilot 放行：** 用户另行批准候选结构生成及六个具体 DFT 标签；D/E/F 环境和运行状态核实、ρref 与 sampler 许可核实、相同 snapshot 的 k 网格/设置冻结、角色/去重账本冻结、独立计算空间可用，并确认单点预算。详细硬门见 `pilot_release_gates.md`。
4. **C 质量审计：** 仅在 owner 归档后核验源/输入/结果 SHA、元素和顺序、PBC、energy/force 单位与键、SCF、有限能量/力、method/mesh、重复、实际 wall/peak-memory、液态多样性。任何 SCF/hash/许可证/数值口径失败均隔离该标签，不进入候选数据表。
5. **正式剩余批次：** 只有用户单独批准后才能处理其余 24 个单点。D/E/F/G/H 不从本文件推导出自动任务授权；本轮停止于计划，不分配运行 owner、不发起作业。

## 依据与复现

上游已发布的数据规划基准为 `3ef38c8af17ed975ea2ceb96baeadce804320eb1`。本报告以当前同步的 `origin/main` `bcb51e5e05493e846ef3bf2cc86081dba8b11abe` 及其 `unified_potential_independent_review_20261010/REPORT.md`、液体 pilot CSV、批次预算 CSV 为依据。审计中，纯相 25 条归档对应 11 个输入几何源；旧角色仍不变。液态 Ag 构型/轨迹未找到，液体密度的 1250 K 精确表值/方程未核实。具体版本/来源记录见 `PROVENANCE.md`，文件校验见 `SHA256SUMS.txt`。

## 独立交接文件状态

为避免各 owner 依赖其他云端的消息，目录另含 `HANDOFF_D.md` 至 `HANDOFF_H.md` 五份完整独立指令、`MASTER_TASK_MATRIX.csv`、`RESOURCE_AND_TIME_BUDGET.md`、`PILOT_RELEASE_GATES.md`、`DATA_IDENTITY_AND_SPLIT.md` 与 `README.md`。所有 30 行都标记 `STRUCTURE_NOT_READY`/无 input SHA；交接文件完整可供只读环境预检，但任何 DFT 仍未获授权。`HANDOFF_READY = YES` 只指文件完整，`DFT_EXECUTION_READY = NO`。
