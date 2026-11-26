# V14 Ag–W 高斯 TTM 加热探索数据

此目录保存 V14 更细双连续 Ag/W 结构及独立的高斯热源试跑。V5–V13 和原始几何数据保留不变。

## 当前结论

- `source_scale=1`（名义 I₀=5 eV/(ps·Å²)）在 2.25 ps 时受热区 Ag 约 283 K，未熔化；为节省计算，将未完成的冷却尾段停止，进度日志保留。
- `source_scale=80`（有效 I₀=400 eV/(ps·Å²)）运行 7.5 ps，软件能量账通过；但 5 ps 关源后仍处于高温，Q6+PTM+MSD 候选在横向近乎均匀分布，没有可识别的局域熔池边界，不能给出熔池尺寸或宣称凝固。
- 初态 Q6/PTM 已高度无序。结构候选数不是液相原子数。完整分析见 [V14_SCALE80_TTM_HEATING_REPORT.md](V14_SCALE80_TTM_HEATING_REPORT.md)。

## Scale=80 主要产物

- `v14_ttm_scale80_fwhm0p5ps_cooling5ps.json`：完整 7.5 ps TTM-MD 能量与温度账。
- `v14_ttm_scale80_fwhm0p5ps_cooling5ps.lammps.log`：LAMMPS 日志。
- `v14_ttm_snapshots_scale80/`：0、0.5、1、2、3、5、7.5 ps 原子快照（GitHub 保存压缩版）。
- `v14_melt_pool_scale80_*.json`：按初态比较的 Q6/PTM/位移候选和径向/深度分布。
- `v14_scale80_t7p5_agw_full_model_disorder_candidates.png`：OVITO 全模型三维图，橙色只表示无序候选，不等于已确认液相。
- `render_v14_scale80_ovito.py` 与 `analyze_v14_melt_pool_q6_ptm_msd.py`：图像与结构分析脚本。

## 历史控制组

- `source_scale=100` 误设高剂量组见 `V14_HEATING_EXPLORATORY_REPORT.md`。
- `source_scale=1` 部分运行的原始进度见 `v14_ttm_scale1_fwhm0p5ps_cooling5ps.progress.json`，该组没有完整结果 JSON。

## 数据约定

坐标单位 Å，时间 ps，LAMMPS metal 单位制。压缩的 dump 快照用 `gzip -dk <file>.gz` 解压。`SHA256SUMS.txt` 为校验清单。
