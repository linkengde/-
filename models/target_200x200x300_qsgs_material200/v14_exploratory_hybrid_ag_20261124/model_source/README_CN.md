# 原始 200 × 200 × 300 Å 双连续模型来源

此模型为 458,248 原子，Ag/Ti₃SiC₂ 质量比约 1:1。`Ag_Ti3SiC2_bicontinuous_200x200x300.data.gz` 是 LAMMPS data 的 gzip 归档；用于 `cases/*` 输入时，需在此目录解压出 `.data` 文件（保留 `.gz`）：

```bash
gzip -k -d Ag_Ti3SiC2_bicontinuous_200x200x300.data.gz
```

`build_model.py`、`qsgs_bicontinuous_mask.npz`、`geometry_report.json` 和 `build.log` 用于追溯结构构建。`extxyz.gz` 是同一结构的 ASE 坐标归档。
