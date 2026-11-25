# V14 Ag–W 高剂量加热探索记录

此目录记录 V14 更细双连续骨架上的一次高剂量 Ag-only Gaussian TTM-MD 试跑。原始 V14/V5–V13 数据均保留。

## 结论

- 7.5 ps 计算及能量账完成，软件耦合账差 1.736 eV。
- 这次误用 `source_scale=100`（有效 I₀=500 eV/(ps·Å²)，为名义 I₀=5 的 100 倍），累计向模型输入约 449.3 keV。
- 末态仍处于高温，没有冷却到固态；Q6/PTM 在初态已显示高度无序，所以不能从这次结果定义可靠的局域熔池尺寸。
- 该结果只作为高剂量探索对照。详见 [V14_HEATING_EXPLORATORY_REPORT.md](V14_HEATING_EXPLORATORY_REPORT.md)。

## 主要文件

- `v14_efs_final.data.gz`：V14 弛豫结构输入。
- `v14_restart_300K.dump.gz`：加热程序读取的初始位置和速度。
- `v14_snapshot_start.dump.gz`、`v14_snapshot_3ps.dump.gz`、`v14_snapshot_7p5ps.dump.gz`：TTM 过程快照。
- `v14_ttm_scale100_fwhm0p5ps_cooling5ps.json`：运行配置、逐时程和能量账。
- `v14_ttm_scale100_fwhm0p5ps_cooling5ps.lammps.log`：LAMMPS 运行日志。
- `v14_melt_pool_q6_ptm_msd_*.json`：初态、3 ps 和 7.5 ps 的结构判据。
- `v14_t7p5_agw_full_model_disorder_candidates.png`：OVITO 全盒三维图；橙色是无序候选，不代表已确认液相。
- `analyze_v14_melt_pool_q6_ptm_msd.py`、`render_v14_disorder_ovito.py`：可复现分析和渲染脚本。
- `SHA256SUMS.txt`：提交数据完整性校验。

压缩轨迹可用 `gzip -dk <file>.gz` 解压。主模型坐标使用 Å、时间使用 ps、LAMMPS 单位制为 metal。
