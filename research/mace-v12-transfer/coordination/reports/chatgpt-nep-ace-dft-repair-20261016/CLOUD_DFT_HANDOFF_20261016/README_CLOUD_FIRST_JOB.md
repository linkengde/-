# 云端DFT交接包 — Ag–Ti₃SiC₂ (2026-10-16)

**已提交到GitHub的真实输入资料**。由**本次ChatGPT聊天界面**整理上传，原DFT标签来自历史GitHub归档，当前没有运行任何新DFT。

## 立即要算的任务（P0）

先解决Ag10Ti12Si2C8、32原子、3D周期的\`Ti_gap_2p3\` k点标签冲突。历史数值：
- 历史1×4×2：native E = −207.463968321659 eV；
- 历史Γ点：native E = −212.001510218357 eV；
- 能量相差4.537541896698 eV，部分原子力分量相差约5.04586 eV/Å。

**两者不同k网格；历史1×4×2是否Γ居中尚未核验，不能直接称为同方法收敛对照。** 
需用同一GPAW/PBE/PW500/PAW/Fermi0.1设置先跑一个新的\`1×4×2\`作为一致网格起点，再依次选择\`1×6×4\`，必要时\`1×8×6\`；两两比对native能、force-consistent自由能和全原子力。禁止一次自动提交三个作业。

### 云端计算所需软件、算力

- Python与已有环境中一致的GPAW **26.7.0**、ASE **3.29.x**、gpaw_data/PAW setup、NumPy、4个MPI ranks，至少对应4个CPU核心；CPU内存/磁盘按实际预检分配，不假设低配够用。
- 设置锁定：**PBE；PW截断500 eV；Fermi-Dirac宽度0.10 eV；不移动任何原子；三维周期边界；自旋极化关闭**。不得私自修改电子数、赝势、晶胞和原子顺序。
- 先核对旧k点中心偏移、原始输入、实际PAW SHA，并保留新的k网格居中方式一致。现有runner显式\`gamma=True\`，与历史的kmesh定义是否一致必须标注。
- 首个作业要求完整的单点能量和逐原子力、SCF信息、GPAW日志、软件版本、PAW SHA和输入SHA。
- 每个网格单独目录。远程云端跑完后先核对结果，不能因为\`SCF converged\`就声明k点收敛。

### 精确执行步骤

\`\`\`bash
git clone --branch audit-nep-ace-dft-reference-repair-20261016 --single-branch https://github.com/linkengde/-.git
cd -
ROOT=research/mace-v12-transfer/coordination/reports/chatgpt-nep-ace-dft-repair-20261016

python "$ROOT/CLOUD_DFT_HANDOFF_20261016/verify_cloud_inputs.py"
python "$ROOT/run_32atom_kmesh_MPIOPT.py" --kmesh 1x4x2

# 以下只在已安装并核实GPAW/ASE/MPI的DFT节点运行：
python "$ROOT/CLOUD_DFT_HANDOFF_20261016/prepare_paw_manifest.py" \
  --output /ABS/OUTSIDE_REPO/paw_manifest.json

# 只启动第一份，并在结果完成后停止：
mpirun -np 4 python "$ROOT/run_32atom_kmesh_MPIOPT.py" \
  --kmesh 1x4x2 --execute \
  --paw-manifest /ABS/OUTSIDE_REPO/paw_manifest.json \
  --output /ABS/OUTSIDE_REPO/Ti_gap_32_k1x4x2
\`\`\`

\`/ABS/OUTSIDE_REPO/\`为云端计算节点实际可写路径的**占位符**，必须替换，绝不当作已建立的目录。不要把\`state.gpw\`、PAW文件、敏感路径或大体积计算缓存直接push进仓库。

完成一份作业后，返回日志、RESULT.json、实际墙钟/内存、收敛状态、SHA，再决定是否继续1×6×4。必须保存原始输出目录，不改写历史DFT。

### P1/P2 已备份的后续数据

\`INPUT_MANIFEST.json\`包含**11个已在GitHub仓库中的无标签extxyz**：1个32原子冲突修复参考；6个135原子Si端Ag/TSC界面（R0/R1，各2.35/2.60/2.95 Å）；2个Ag32 FCC静态扰动；2个Ti₃SiC₂ 48原子静态扰动。其SHA是**Git blob SHA-1**，不是DFT总能和不是SHA-256；预检会同时计算真正的文件SHA-256。

任务优先级在\`DFT_JOB_QUEUE.csv\`。P1/P2的kmesh/真空/dipole/spin尚需分别核对，**禁止直接套32原子的1×4×2/1×6×4设置运行135原子slab**。

这些输入不能解决液态Ag、高温MAX/液Ag–MAX和固液共存的缺失；Ag64非DFT轨迹只有在其真正有证据证明处于液态、映射完整、完成DFT单点时才可用于训练。当前无该轨迹在本交接包内的文件。

### 数值收敛判据（事先冻结）

1. 相邻相同物理设置的kpoint网格之间native总能差 / 32 ≤ **2 meV/atom**；
2. force-consistent自由能差 / 32 ≤ **2 meV/atom**；
3. all-atom力向量RMS差 ≤ **0.01 eV/Å**，同时报告最大逐原子力矢量差；
4. 不得将历史不同method设置/不同spin/不同PAW结果混作收敛证明；
5. DFT结果包含带索引原子顺序、完整forces(N×3)和输入哈希。SCF未收敛标记\`INVALID_LABEL\`，不用于NEP/ACE训练。

### 必须返回

\`DFT_P0_RESULT.json\`（或runner产生的\`RESULT.json\`）、\`gpaw.log\`、\`SCF_AND_KPOINT_AUDIT.md\`、\`ENVIRONMENT_AND_PAW_SHA.json\`和输入/结果文件SHA-256清单；存独立目录，push新Git提交并fresh-fetch核对。若缺权限、依赖、PAW或者内存，不要降参数凑结果，返回\`BLOCKED\`和具体失败依据。

**源头**：GitHub \`linkengde/-\`，分支 \`audit-nep-ace-dft-reference-repair-20261016\`。本目录是ChatGPT聊天直接生成的交接资料，并非已生成的新DFT结果。
