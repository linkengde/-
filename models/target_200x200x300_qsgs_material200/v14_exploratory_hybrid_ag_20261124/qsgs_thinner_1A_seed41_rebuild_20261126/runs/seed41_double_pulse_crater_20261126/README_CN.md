# seed41 第二、三次局域脉冲归档

本目录补充第一脉冲后 5.8–14.8 ps 的延长计算。完整结论、边界条件和文件说明见 [seed41_three_pulse_report_CN.md](seed41_three_pulse_report_CN.md)。第一脉冲阶段记录在相邻目录 `../seed41_local_pulse_crater_20261126/`。

## 主要图片

- [14.8 ps 全模型 OVITO 三维图](seed41_three_pulse_final_14p8ps_OVITO_CN_300dpi.png)
- [14.8 ps 中心区域 OVITO 图](seed41_three_pulse_crater_closeup_14p8ps_OVITO_CN_300dpi.png)
- [三脉冲形貌、Ag 高空计数、MSD 和 Q6 综合图](seed41_three_pulse_diagnostics_CN_300dpi.png)
- [第二脉冲阶段诊断图](seed41_double_pulse_diagnostics_CN_300dpi.png)

## 复现顺序

在本目录执行。归档中的 `.gz` 文件是为仓库体积压缩的原始文本轨迹、dump 或 LAMMPS restart；先恢复所需文件名：

```bash
gzip -dk seed41_pulse2_1479eV_start.dump.gz
gzip -dk seed41_after_pulse2_7p8ps.restart.gz
gzip -dk seed41_afterpulse2_10p8ps_nomsd.restart.gz
gzip -dk seed41_pulse3_1479eV_start.dump.gz
gzip -dk seed41_after_pulse3_12p8ps.restart.gz
```

Pulse2 初始 restart 与 5.8 ps 输入 dump 位于相邻的 `seed41_local_pulse_crater_20261126` 目录。

```bash
python prepare_seed41_pulse2.py
lmp -sf omp -pk omp 4 -in in.seed41_pulse2_1479eV_2ps.lmp
lmp -sf omp -pk omp 4 -in in.seed41_afterpulse2_nve_7p8_to10p8ps_nomsd.lmp
python prepare_seed41_pulse3.py
lmp -sf omp -pk omp 4 -in in.seed41_pulse3_1479eV_2ps.lmp
lmp -sf omp -pk omp 4 -in in.seed41_afterpulse3_nve_12p8_to14p8ps_nomsd.lmp
```

两段较长 NVE Ag 轨迹拆成 `part01` 等独立 gzip 分片以满足仓库单文件限制；若要恢复完整轨迹，按文件名顺序逐片解压并串接，例如：`for p in seed41_afterpulse2_NVE_ag_nomsd.part*.dump.gz; do gzip -cd "$p"; done > seed41_afterpulse2_NVE_ag_nomsd.dump`。对 pulse3 后的 NVE 轨迹使用对应 `seed41_afterpulse3_NVE_ag_nomsd.part*.dump.gz`。

绘图与后处理还需要 Python `numpy`、`matplotlib`、Pillow 和 OVITO Python 模块。图片中的红色表示 Ag 位于原始材料表面 z>220 Å，不等于已蒸发；橙色表示 PTM 未分类为 FCC/HCP，也不单独等于熔化。
