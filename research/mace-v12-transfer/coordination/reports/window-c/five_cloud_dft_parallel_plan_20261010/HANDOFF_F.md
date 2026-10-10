# 云端 F：Ag–Ti₃SiC₂ 五云端 DFT 交接

项目：Ag–Ti₃SiC₂ 双连续复合触头统一 Ag–Ti–Si–C MACE 势开发。最终应用包含电弧熔池、Ag/骨架高温作用、熔池冷却凝固和后续熔坑/喷溅研究。

**状态：交接材料可供主机预检；DFT 执行未放行。** 收到这份文件不构成结构生成、DFT、MACE 训练、MD/AIMD、LAMMPS 或 TTM 授权。本轮只按这里列出的六个 task ID 工作，不接手其他云端任务，不修改 A/B。

依据仓库 `linkengde/-` 的已归档规划提交 `3ef38c8af17ed975ea2ceb96baeadce804320eb1` 和完成状态提交 `81e7a218312f187957538cda322bf49b7e9dceae`。开始后先读取当前 `origin/main` 中 `research/mace-v12-transfer/coordination/reports/window-c/unified_potential_dft_data_plan_20261010/` 下的 `REPORT.md`、`proposed_batches.csv`、`METHOD_AND_ROLE_RULES.md`、`split_plan.csv`、`existing_data_inventory.csv` 和 `PROVENANCE.json`；确认两个指定提交在同步历史中，不 reset 到旧快照。

## 本云端六项任务

所属云端：**F**。本云端唯一获分配的提案 ID：`C-F-AGL-1050K-rho103-seed01`, `C-F-AGL-1250K-rho100-seed02`, `C-F-AGL-1250K-rho103-seed01`, `C-F-AGL-1600K-rho097-seed01`, `C-F-TSC-PF01-CFG03`, `C-F-AG-PF01-CFG03`。不得启动或重排其他云端的 task ID。全部行当前为 `STRUCTURE_NOT_READY`：没有可运行坐标、DFT 输入或输入 SHA。parent/seed 数值均为计划值，source hash 与实际 RNG 轨迹尚未产生。

| Task ID | 系统与原子数 | 温度/密度或构型 | Proposed parent / seed | 来源、输入和 SHA 状态 | 预计成本 |
|---|---|---|---|---|---:|
| `C-F-AGL-1050K-rho103-seed01` | liquid_Ag_candidate；64 atoms | 1050 K; 1.03 × rho_liquid(T); 绝对密度来源待核; independent trajectory/seed snapshot after liquid-state, RDF, coordination, order, MSD and diversity checks | `PF-AGL-1050K-rho103-seed01`; seed `105010301` (proposed) | Future geometry-only sampling; approved sampler not selected. Candidate MACE/EAM energies and forces are not labels.; NOT_GENERATED; source, sampler/license and density basis pending; STRUCTURE_NOT_READY; no runnable input exists; SHA: NOT_COMPUTED; generate only after separate structure-generation approval | 3-8 h / 4 MPI ranks |
| `C-F-AGL-1250K-rho100-seed02` | liquid_Ag_candidate；64 atoms | 1250 K; 1.00 × rho_liquid(T); 绝对密度来源待核; independent trajectory/seed snapshot after liquid-state, RDF, coordination, order, MSD and diversity checks | `PF-AGL-1250K-rho100-seed02`; seed `125010002` (proposed) | Future geometry-only sampling; approved sampler not selected. Candidate MACE/EAM energies and forces are not labels.; NOT_GENERATED; source, sampler/license and density basis pending; STRUCTURE_NOT_READY; no runnable input exists; SHA: NOT_COMPUTED; generate only after separate structure-generation approval | 3-8 h / 4 MPI ranks |
| `C-F-AGL-1250K-rho103-seed01` | liquid_Ag_candidate；64 atoms | 1250 K; 1.03 × rho_liquid(T); 绝对密度来源待核; independent trajectory/seed snapshot after liquid-state, RDF, coordination, order, MSD and diversity checks | `PF-AGL-1250K-rho103-seed01`; seed `125010301` (proposed) | Future geometry-only sampling; approved sampler not selected. Candidate MACE/EAM energies and forces are not labels.; NOT_GENERATED; source, sampler/license and density basis pending; STRUCTURE_NOT_READY; no runnable input exists; SHA: NOT_COMPUTED; generate only after separate structure-generation approval | 3-8 h / 4 MPI ranks |
| `C-F-AGL-1600K-rho097-seed01` | liquid_Ag_candidate；64 atoms | 1600 K; 0.97 × rho_liquid(T); 绝对密度来源待核; independent trajectory/seed snapshot after liquid-state, RDF, coordination, order, MSD and diversity checks | `PF-AGL-1600K-rho097-seed01`; seed `160009701` (proposed) | Future geometry-only sampling; approved sampler not selected. Candidate MACE/EAM energies and forces are not labels.; NOT_GENERATED; source, sampler/license and density basis pending; STRUCTURE_NOT_READY; no runnable input exists; SHA: NOT_COMPUTED; generate only after separate structure-generation approval | 3-8 h / 4 MPI ranks |
| `C-F-TSC-PF01-CFG03` | Ti3SiC2_skeleton；48-96 atoms | static / controlled small displacement; source ionic temperature not selected; 固态；cell density 待 source cell 冻结; C_mode_displacement | `PF-TSC-PF01`; seed `7310103` (proposed) | New supercell/orientation or independent thermal/mode parent from a provenance-reviewed structure; no geometry generated; NEW_PARENT_NOT_SELECTED; exact coordinates/hash pending; exclude scored parent families; STRUCTURE_NOT_READY; no runnable input exists; SHA: NOT_COMPUTED; generate only after separate structure-generation approval | 3-6 h / 4 MPI ranks |
| `C-F-AG-PF01-CFG03` | solid_Ag；32-64 atoms | static / controlled small response; source ionic temperature not selected; 固态；cell density 待 source cell 冻结; small_shear_or_mode | `PF-AG-PF01`; seed `7410103` (proposed) | New fcc Ag supercell/seed from a provenance-reviewed parent, structurally independent of scored numerical controls; NEW_PARENT_NOT_SELECTED; exact coordinates/hash pending; exclude scored parent families; STRUCTURE_NOT_READY; no runnable input exists; SHA: NOT_COMPUTED; generate only after separate structure-generation approval | 2-4 h / 4 MPI ranks |

角色：全部是未来训练候选的提案，逐条仍须数据来源、energy/force 和角色审核；不是独立开发集或测试集。原归档中的 15 个纯 Ag、10 个 Ti₃SiC₂ 数值控制标签继续保持原角色。

## pilot 阶段状态

本云端分到六个液体 pilot 的 2 个候选：`C-F-AGL-1250K-rho100-seed02`, `C-F-AGL-1250K-rho103-seed01`。它们就是上表中 30 个任务的子集，**不是新增结构**。两项运行也必须等待用户对真实结构生成和 DFT 的明确批准。

1250 K Ag 熔点附近参考液态密度的精确文献表/方程尚未核实，六个 64-atom 液体坐标/轨迹也尚未生成。液体采样器及许可尚未冻结。不得把计划参数补写成已有坐标或标签。结构生成如需 MD/AIMD，须另有明确授权。

## 第一阶段：只读预检查（不执行 DFT）

1. **资源：** 报告物理/逻辑 CPU 核、当前可用内存、计算工作盘剩余、Linux/系统版本、当前 MPI/GPAW/ASE/Python/PAW 数据版本和许可状态、BLAS/LAPACK/ScaLAPACK/LibXC 构建信息。目标参考：一个 4-rank 作业至少 4 个可用物理核；64-atom pilot 至少 16 GiB 可用内存；较大正式 cell 建议 32 GiB；每个独立工作区建议 ≥30 GB。它们是筛选建议，不是已经测出的主机性能。
2. **可复制的只读采集命令：**

```bash
python3 - <<'PY'
import os, sys, shutil, importlib.metadata as md
from pathlib import Path
print('python', sys.version.split()[0], 'logical_cpu', os.cpu_count())
try:
    blocks=Path('/proc/cpuinfo').read_text().split('\n\n')
    pairs=set()
    for block in blocks:
        d=dict(line.split(':',1) for line in block.splitlines() if ':' in line)
        if 'physical id' in d and 'core id' in d: pairs.add((d['physical id'].strip(),d['core id'].strip()))
    print('physical_cpu', len(pairs) or 'unknown')
except Exception: print('physical_cpu unknown')
try: print('memory', Path('/proc/meminfo').read_text().split('MemAvailable:')[1].splitlines()[0].strip())
except Exception: print('memory unknown')
st=os.statvfs('.'); print('workdir_free_GB', round(st.f_bavail*st.f_frsize/1e9,1))
for p in ('gpaw','gpaw-data','ase','numpy','scipy'):
    try: print(p, md.version(p))
    except md.PackageNotFoundError: print(p, 'not-installed')
print('mpi_launcher', bool(shutil.which('mpiexec') or shutil.which('mpirun')))
print("GPAW_SETUP_PATH configured", bool(os.environ.get("GPAW_SETUP_PATH")))
try:
    import gpaw_data
    print('PAW data directory exists', Path(gpaw_data.datapath()).is_dir())
except Exception: print("PAW data import/directory unavailable")
PY
if command -v mpiexec >/dev/null 2>&1; then mpiexec --version | sed -n '1,2p'; fi
uname -srm
git fetch origin main
git status --short --branch
git rev-parse origin/main
git merge-base --is-ancestor 3ef38c8af17ed975ea2ceb96baeadce804320eb1 origin/main
git merge-base --is-ancestor 81e7a218312f187957538cda322bf49b7e9dceae origin/main
```

不要粘贴凭证、`env` 全量输出或私人绝对路径。MPI binary 存在不等于与 GPAW 的 MPI ABI 兼容；由本机环境 owner 核对 GPAW MPI build 与 runtime、线性代数库。核实 PAW data 1.2.1 的来源、许可、可用性及文件 SHA，不要擅自安装/更换组件或科学参数。

3. **仓库/任务：** 记录当前 main SHA 和本地 branch/status；只从当前已发布 main 安全快进。使用云端自己的 branch（建议 `cloud-f-b1-dft`）和自己的 `${DFT_WORK_ROOT}`。先只读报告已有、运行中、排队、未启动 task；活动计算和 checkpoint 完全保持原样。输入列表目前都不存在，应记录 `STRUCTURE_NOT_READY`，不能生成 DFT 输入。
4. **结构/输入：** 未获得 C 冻结的结构包前，不创建 task 目录/坐标/DFT 输入。日后拿到输入后，先核 task ID 唯一性、来源记录、parent/seed、原子数、元素与 atom ID/排列、晶胞、PBC、温度/密度、几何哈希和 input SHA；与所有已公开非封存结构和其他云任务做精确及周期近重复核对。

## 科学设置（仅用于后续获批作业）

- GPAW 26.7.0 + ASE 3.29.0 + PBE + PAW data 1.2.1；PW cutoff 500 eV；Fermi 电子展宽 σ=0.10 eV；SCF density/energy/eigenstate 阈值沿用经审查的 1e-5 输入语义。不得自行更改。
- 保存 native `energy`、与计算力一致的 `free_energy`，以及每个原子的三维力 `REF_forces`。能量单位 eV、力 eV/Å；训练候选能量键和力口径必须配对。另存 cell/PBC、k 网格、计算版本及所有参数。
- **液态 Ag：** 同一个冻结 64-atom snapshot 对照 Γ 与 2×2×2；若 ΔE>2 meV/atom、force RMS>0.01 eV/Å 或最大力差>0.03 eV/Å，扩到 3×3×3 或更密并报 C。未满足阈值不放行本状态点的正式标签；绝不全批默认 Gamma。
- **Ti₃SiC₂：** 按原始已审计骨架网格的倒格矢密度缩放到最终新超胞（primitive reference 包含 10×10×4）；至少对一个代表性新 parent 比较匹配密度网格与更密网格，按上列差异门验收。晶胞未冻结则不能给数值 k 网格。
- **固态 Ag：** 依据已审计 Ag parent 的收敛倒格矢密度映射到新超胞；同一代表性结构比较目标密度与更密网格，按同一差异门验收。Ag 为金属，不能未经测试一律 Gamma。
- mesh 预算只是预注册目标，不是精度结果。方法、PBC、密度、元素顺序、参数口径或 mesh 不匹配时，停下汇报，不擅自改参数。

## 分阶段执行与验收

1. **现在：** 仅完成只读环境、任务状态、数据来源、许可和输入存在性报告。不做结构生成、k 网格 DFT pilot 或正式 DFT。
2. **用户明确批准 pilot 后：** 仅执行本文件标成首批 pilot 的任务：`C-F-AGL-1250K-rho100-seed02`, `C-F-AGL-1250K-rho103-seed01`。pilot 候选仍须等真实液体结构生成及 density/sampler 来源解锁。
3. **C 审查 pilot 后：** 提交每个结果的 SCF 收敛、energy/free-energy、完整原子力、k 网格对照、结构多样性/液态判据、SHA、verifier、wall time 和 peak memory。等待 C 明确接受或退回。
4. **剩余正式单点：** 只有用户进一步批准，才可执行本云端上表剩余 HOLD 项；未经批准不以 pilot 通过为自动放行。

每个获批作业暂定单例上限 12 h wall 或 160 SCF iteration；首例超过 8 h 先停批报 C 复估。保存原日志和 checkpoint，不自动重启。验收必须 SCF `converged=true`、能量/力全部有限、atom order/PBC/hash 正确、mesh 差满足预算且 archive verifier PASS。

## 独立目录、归档和停止条件

作业获批后，每个 task 在本云端独立路径 `${DFT_WORK_ROOT}/f/<task_id>/{input,checkpoint,logs,output}`。每台云端的路径和 checkpoint 不共享。至少归档冻结输入结构及 source manifest、DFT 参数、完整 SCF log、native energy、free energy、原子力、cell/PBC/order、硬件与运行耗时/内存、input/result SHA256、verifier 状态/失败理由。不要提交活跃 checkpoint 或 force-push main。

出现结构不存在/来源不明、重复、原子映射或 PBC 错误、PAW/许可不明、资源不足、超预算、checkpoint 不一致、NaN/Inf、SCF/mesh 未达标或参数不符时：停止该 task，保留日志和 checkpoint，报告 C；不静默改几何、k 点、展宽或 cutoff 后继续。

## 本云端交接状态

- `HANDOFF_READY = YES`：表示文件和任务清单完整，可以执行上述只读预检。
- `DFT_EXECUTION_READY = NO`：上表所有结构/输入/hash 未就绪；用户尚未放行任何计算。
- 当前 owner 回传脱敏的资源、环境、任务运行状态与来源缺口后，等待 C 复核和用户明确批准。
