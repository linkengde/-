# Round7 347,464 原子 TTM 加热与熔化判读复现包

本目录保存一次完整的 **3 ps 探索性热源开启试跑**、轨迹帧、能量账本，以及用 Q₆、PTM、MSD 联合判读的脚本和报告。3 ps 结果显示热区明显失序，但不足以确认熔化。随后同条件长时试跑在 5.78 ps 停止；局部熔化判读和扩展轨迹见相邻目录 `../round7_347k_melt_onset_5p78ps/`。

长时试跑配置与本次 3 ps 配置完全相同，只把 `AG_TTM_TRIAL_STEPS` 从 6,000 改为 16,000。它从原始几何和相同随机种子重新开始，并非从 3 ps 末态接续：现有 driver 没有保存 TTM 电子网格检查点，不能无损恢复耦合状态。长时试跑在 5.78 ps 达到停止判据后主动结束；5.5 ps 是最后一个完整轨迹帧。

## 复现输入和来源

- 模型：347,464 原子；Ag 80,674，Ti₃SiC₂ 266,790；盒子 200 × 200 × 300 Å，材料占 z=0–150 Å。
- 归档来源分支：`codex/round7-v10-artifacts-20261118`，提交 `36b203301e05b981139016d540aa1b36d6929c93`。
- 原始 LAMMPS data 文件 SHA-256：`76b3ec5b2ee0a233dc7f69d760c0d7f5733bd899e8e833f60121553c8296606f`。
- Ag–Ag：Foiles `Ag_u3.eam`，文件 SHA-256：`cb3d605cd27ec6efb50c4dd7462c2ab3d967e63775dc4400de9134315a44d807`。
- Ag–骨架：WCA 排斥几何边界，参数见 `potential.geometry_wca_fixed_skeleton.inc`。
- TTM：仓库中的 Python 参考算子 `ag_ttm_reference.py` 经 LAMMPS Python API 耦合；不是 LAMMPS 原生 `fix ttm`。
- 温度边界：所有 Ti/Si/C 原子固定，Ag 初始 z<20 Å 固定；本加热段没有单独恒温热沉。
- 初速度：移动 Ag 以 300 K、随机种子 190427 初始化。
- 数值条件：0.5 fs，4 个 OpenMP 线程，6,000 步（3 ps）；XY 圆对称高斯 × z 向表面深度指数衰减，R=30 Å、skin=20 Å，峰值 5 eV/(ps·Å²)。

## 获取结构并核验

在仓库根目录 `/workspace/-` 执行。结构本体保留在上述归档分支中，避免在本结果目录重复存储 28 MB 的原始 data 文件。

```bash
git fetch origin codex/round7-v10-artifacts-20261118
git show 36b203301e05b981139016d540aa1b36d6929c93:archive/round7-v10-artifacts-20261118/round7_pod_snap/models/target_200x200x300_qsgs/Ag_Ti3SiC2_bicontinuous_200x200x300.data > /tmp/Ag_Ti3SiC2_bicontinuous_200x200x300.data
echo '76b3ec5b2ee0a233dc7f69d760c0d7f5733bd899e8e833f60121553c8296606f  /tmp/Ag_Ti3SiC2_bicontinuous_200x200x300.data' | sha256sum -c -
```

## 软件环境

原试跑环境使用 Python 3.12、LAMMPS Python 包 2025.7.22.4.0、mpi4py 4.1.2、NumPy 2.5.3 和 SciPy 1.18.1。轨迹分析使用 OVITO Python 3.17.2。LAMMPS、TTM 参考代码、EAM 和 WCA 输入分别位于仓库的 `models/ag_ti3sic2_qsgs_1to1_finer/` 下。

## 运行命令

先为新运行目录复制几何元数据，避免覆盖本目录的归档结果。下方命令假设云端虚拟环境位于 `/workspace/.venvs/qsgs`，仓库位于 `/workspace/-`。

```bash
repo=/workspace/-
model=/tmp/Ag_Ti3SiC2_bicontinuous_200x200x300.data
out=/workspace/round7_347k_replay_3ps
mkdir -p "$out"
cp "$repo/models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/geometry_constraint_B/ttm_trial_20261125/round7_347k_3ps/geometry_report.json" "$out/"
cp "$repo/models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/geometry_constraint_B/ttm_trial_20261125/round7_347k_3ps/ttm_round7_347k.py" "$out/"
source /workspace/.venvs/qsgs/bin/activate
export LD_LIBRARY_PATH=/workspace/.venvs/qsgs/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=4
export AG_TTM_TRIAL_STEPS=6000
export AG_TTM_SOURCE_SCALE=1.0
export AG_TTM_TRIAL_OUTPUT_DIR="$out"
export AG_TTM_MODEL_DATA="$model"
python -m py_compile "$out/ttm_round7_347k.py"
python -u "$out/ttm_round7_347k.py" > "$out/run.log" 2>&1
```

要重跑长时长方案，将 `AG_TTM_TRIAL_STEPS=6000` 改成 `AG_TTM_TRIAL_STEPS=16000`，并给 `AG_TTM_TRIAL_OUTPUT_DIR` 指定另一个空目录。原运行在 5.78 ps 根据 Q₆/PTM/MSD 的持续液态样信号主动停止；其他热源参数不变。

若仓库不在 `/workspace/-`，请同步调整脚本内 `ROOT` 常量，使其指向 `models/ag_ti3sic2_qsgs_1to1_finer`。

## Q₆、PTM、MSD 分析

将本目录的 `analyze_q6_ptm_msd.py` 与 `geometry_report.json` 复制到实际输出目录，再使用 OVITO Python 环境执行：

```bash
cp /workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/geometry_constraint_B/ttm_trial_20261125/round7_347k_3ps/analyze_q6_ptm_msd.py "$out/"
source /workspace/.venvs/ovito/bin/activate
python "$out/analyze_q6_ptm_msd.py"
```

分析固定跟踪初始 r_xy≤30 Å 且 z≥120 Å 的 1,078 个 Ag 原子；Q₆/q̄₆ 使用 Ag-only 邻居和 3.6 Å cutoff；PTM 只对 Ag 子集运行；MSD 用初始帧为原点，并报告去除热区整体平移后的值。q̄₆=0.45 仅作同模型趋势筛查，不是普适固液阈值。

## 本次运行结论

3 ps 结束时，全运动 Ag 平均温度 803.75 K，热源加权晶格温度 780.87 K。热区 q̄₆ 中位数由 0.578 降到 0.280，PTM-FCC 比例由 7.70% 降至 0.93%，去整体漂移 MSD 增至 4.739 Å²；这些指标说明明显失序和迁移，但初始结构未弛豫且初始 PTM-FCC 识别率低，温度也未进入块体 Ag 熔化区间，因此不能据此宣称形成了可靠液态 Ag。

本模型的 Ti₃SiC₂ 是完全固定的刚性几何支撑，不振动、不导热；WCA 也不是已拟合的 Ag/Ti₃SiC₂ 物理界面。结果仅用于检查 TTM 数值耦合与结构判读流程，不能作为完整触头熔化—凝固模拟结果。

各文件 SHA-256 记录在 `SHA256SUMS.txt`。轨迹与完整运行指标在本目录归档；原始结构按上面的归档提交取回。
