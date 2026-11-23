# Round7 论文图件第一阶段：图2与图4

本目录只包含已有 Round7 轨迹的后处理结果，不启动新 MD，不修改原始轨迹、势函数或已有计算结果。PNG/TIFF 为 300 dpi；PDF 用于排版和矢量查看。

## 交付文件

- `Figure2_round7_fiveframe_structure_evolution.png/.tif/.pdf`：五帧完整原子模型加 1.60 ps 热源邻域放大。Ti、Si、C、Ag 分别着色；完整模型使用相同正交相机、盒体尺寸和原子半径，并标注 50 Å / 20 Å 物理标尺。
- `Figure4_round7_ag_disorder_candidate_region.png/.tif/.pdf`：五帧 Q6 均值、PTM 非 FCC 持续比例、重算的 Ag 配位数和最大联合候选连通分量。
- `round7_stage1_quality_report.json`：分支核查、源文件 SHA256、原子 ID/类型对齐、候选计数、连通分量复算、配位数及运行异常。
- `figure_captions_and_metrics.md`：中文图注、物理量定义、单位、数据来源和限制。
- `render_round7_stage1.py`：可重复运行的 OVITO/Matplotlib 后处理脚本，含中文说明和注释。

## 数据来源与分支判断

已检查 `main`、`codex/round7-main-pool-figure-20261121`、`codex/round7-v12m-cold-diagnostics-20261119` 和 `codex/upload-round7-fiveframe-trajectory-1p6ps-20261109`。前两类其他分支分别只有上一版图件、单帧候选图或势函数冷态诊断，没有本次所需的五帧完整轨迹。图2因此使用 11 月 9 日轨迹归档；图4使用同一归档中的**更正径向源区**分析目录，避免旧版筛选遗漏 y 条件的问题。

来源仓库：`linkengde/-`，数据提交：`878328378869ba4e5860a4081420fcdc59966c52`。

图2输入：

`models/ag_ti3sic2_qsgs_200x200x150_md_finer_skeleton/runs/round7_fiveframe_joint_analysis_to_1p6ps_20261109/`

图4输入：

`models/ag_ti3sic2_qsgs_200x200x150_md_finer_skeleton/runs/round7_corrected_radial_core_analysis_1p6ps_20261109/`

当前 `main` 不含这些 43 MB 源轨迹文件；运行脚本前需从上述数据分支检出/读取输入，并传入五帧归档和更正径向目录中的文件。SHA 清单会校验轨迹、运行摘要、CSV、诊断 JSON、候选 dump、最大分量 dump 和筛选 JSON。

## 运行方式

```bash
python render_round7_stage1.py \
  --trajectory /path/to/round7_5frames_full_system_0p6_to_1p6ps.dump.gz \
  --source-core-csv /path/to/round7_1p6ps_source_core_ag.csv \
  --ptm-diagnostic /path/to/ptm_coordination_diagnostic.json \
  --candidates-dump /path/to/round7_1p6ps_joint_candidates.dump \
  --largest-dump /path/to/round7_1p6ps_largest_candidate_cluster.dump \
  --screen-json /path/to/round7_1p6ps_joint_screen.json \
  --run-summary /path/to/run_summary_1p6ps.json \
  --trajectory-manifest /path/to/fiveframe_SHA256SUMS.json \
  --analysis-manifest /path/to/corrected_radial_SHA256SUMS.json \
  --out-dir .
```

脚本依赖 OVITO Python 3.17.2、NumPy、Matplotlib 和 Pillow。环境中的 OVITO StandardRenderer 需要 Vulkan；本次在共享环境用 Mesa lavapipe 运行。

## 结果解释范围

874 个联合筛选原子是 Q6/PTM/净位移阈值下的探索性候选，不等同于已确认液态或熔池。阈值来自 Round5，未针对 Round7 校准。Round7 使用开发版 MTP 和 Ag-only Python 参考 TTM；完整报告记录最大力和能量残差异常。图件不能单独证明熔化、熔焊、骨架破坏或界面热传递机制。
