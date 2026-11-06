# P1P Ag70 短时运行归档

该运行从 300 K 初态起步，原设置的最大时长为 7 ps。为做 P1P 源项短跑检查，进程在 3.0128048117 ps 检查点完成后停止；日志中还记录了检查点后的一个短步至 3.0527 ps。它没有完成 7 ps、目标温度保持或冷却。

## 计算设置

- 热源：Ag–Ti₃SiC₂ 复合热源，4.0 mW。
- Marangoni：切向倍率 3×；原法向表面张力未放大。
- 电弧压力：探索性峰值 501 MPa。
- 蒸发：P1P 开启，Wang et al. (2024) 候选 `p_sat(T)`、由其斜率派生的 `L_v`。
- 局部 Ag 蒸气分压：0 Pa，自由蒸发上限情景，非实测工况。
- α：Table 8 锚点 1373–1773 K 内插值；范围外使用明确的敏感性设定 0.1428。
- 网格：64×64×76，单元尺寸约 0.301857 nm；`dt_max=20 fs`。
- 顶边界开放至环境压力；不解蒸气域、不算凝结或反冲压力。

## 结果文件

- `p1h_history.csv`：温度、液态 Ag 体积、自由表面统计、能量和时间步历史。
- `p1p_evaporation_flux.csv`：界面面积、最大通量、α 来源计数、累计逸出质量和账本残差。
- `p1j_transport_validation.csv`、`p1h_arc_pressure_history.csv`、`p1h_peak_pool_dimensions.csv`：耦合和压力/熔池诊断。
- `p1p_real_7ps.stderr`：完整运行日志；名称保留原运行输出名，实际运行在约 3 ps 停止。
- `run_config.env`：参数摘要。
- `checkpoint/`：压缩的最后完整 Basilisk dump 及其配套状态文件。

## 恢复检查点

在具有相同运行代码、网格与输入几何的工作目录中，解压 dump 并保留同名 sidecars：

```bash
gzip -dc checkpoint/p1h_restart_i0000000161_t003.012805_ps.dump.gz \
  > p1h_restart_i0000000161_t003.012805_ps.dump
cp checkpoint/p1h_restart_i0000000161_t003.012805_ps.dump.state .
cp checkpoint/p1h_restart_i0000000161_t003.012805_ps.dump.p1p .
cp checkpoint/p1h_restart_i0000000161_t003.012805_ps.dump.complete .
```

原始 dump SHA256：`b0059aaf6384435ebea87f8859fe62749772521702a2b310631086b6d06b2e1f`。

先用 `sha256sum -c SHA256SUMS.txt` 验证归档。该重启点保留的是约 3.0128 ps 状态，不代表该算例已冷却或完成。
