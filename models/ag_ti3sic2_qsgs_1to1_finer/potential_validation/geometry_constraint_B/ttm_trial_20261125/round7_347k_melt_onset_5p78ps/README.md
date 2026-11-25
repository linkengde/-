# 347k TTM 长时段：局部 Ag 熔化判读

本目录归档从相同 300 K 初始构型重新开始的长时 TTM 试跑。计划时长为 8 ps（16,000 步），在 5.78 ps 根据热区 Q₆、PTM 与去漂移 MSD 的持续液态样信号主动停止。**实际最后一条积分记录是 5.78 ps；最后一个完整轨迹帧是 5.5 ps。** 轨迹未执行冷却或凝固。

完整前 0–3 ps 轨迹帧位于兄弟目录 `../round7_347k_3ps/`；本目录只存放新增的 3.25–5.5 ps 帧，以减少重复文件。`analyze_q6_ptm_msd.py` 会自动合并两处帧后分析 0–5.5 ps。

## 运行条件

长时试跑和 3 ps 试跑使用同一模型、Ag EAM、WCA 几何边界、随机种子 190427、TTM 参数及 0.5 fs 步长；只有计划步数由 6,000 改为 16,000。具体模型来源、模型 SHA-256、势函数、软件版本和 3 ps 完整复现命令见 `../round7_347k_3ps/README.md` 与 `DEPENDENCY_SHA256.txt`。

运行使用的参数：

```bash
export AG_TTM_TRIAL_STEPS=16000
export AG_TTM_SOURCE_SCALE=1.0
export OMP_NUM_THREADS=4
export AG_TTM_TRIAL_OUTPUT_DIR=/workspace/round7_347k_long_replay
export AG_TTM_MODEL_DATA=/tmp/Ag_Ti3SiC2_bicontinuous_200x200x300.data
```

需要把 `geometry_report.json` 和 `ttm_round7_347k.py` 放入输出目录，再按兄弟 3 ps README 中的命令运行。脚本不会在线计算 Q₆/PTM/MSD；长时方案是先逐 0.25 ps 输出帧，再离线分析，确认三项信号持续后停止热输入。

## 分析复现

使用兄弟目录的初始至 3 ps 帧和本目录 3.25–5.5 ps 帧。分析器固定追踪初始 r_xy≤30 Å、z≥120 Å 的 1,078 个 Ag 原子。Q₆/q̄₆ 是 Ag-only、3.6 Å cutoff；q̄₆>0.45 只作同模型趋势筛查。MSD 使用初始帧并去除热区整体平移。

```bash
source /workspace/.venvs/ovito/bin/activate
python analyze_q6_ptm_msd.py
```

Ag–Ag 最近邻距离另由 `analyze_ag_ag_neighbors.py` 计算，需使用含 SciPy 的模拟 Python 环境：

```bash
source /workspace/.venvs/qsgs/bin/activate
python analyze_ag_ag_neighbors.py
```

结果表为 `q6_ptm_msd_timeseries.csv` 和 `msd_timeseries.csv`，完整逐帧方法说明在兄弟目录中的分析脚本与摘要 JSON。

## 结论与限制

5.0–5.5 ps 期间，热区 q̄₆ 中位数稳定在 0.235，PTM-FCC/HCP 都降至 0，去漂移 MSD 增至 10.168 Å²（RMS 3.189 Å）。这支持热源区 Ag 的局部液态样熔化；不代表整块 Ag 熔化。初始结构未弛豫，初始 PTM-FCC 仅 7.70%，温区均值仍低于块体 Ag 熔点附近，因此这是探索性结构判定，不是平衡熔点或真实触头熔池的定量预测。

5.5 ps 最近 Ag–Ag 邻距最小 2.220 Å，1% 分位数 2.400 Å；虽然没有达到 50 eV/Å 安全力阈值或出现原子丢失，这类短接触仍需保留为风险信号。骨架完全固定，WCA 界面无真实 Ag/MAX 黏附，TTM 使用 Python 参考实现。

最终数据、热量账本及判读见 `MELTING_ASSESSMENT_Q6_PTM_MSD_5P78PS.md`、`TRIAL_SUMMARY_STOPPED_5P78PS.json`、`step_metrics.csv` 和 `run.log`。文件校验和在 `SHA256SUMS.txt`。
