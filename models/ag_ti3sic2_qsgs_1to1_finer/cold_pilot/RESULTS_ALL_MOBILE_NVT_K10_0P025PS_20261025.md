# 全活动区 NVT 的冷态预平衡诊断：0.025 ps

日期：2026-10-25（Asia/Shanghai）。该诊断检查将所有活动原子置于 300 K NVT、并用 K=10 eV/Å² 弱弹簧约束活动骨架，能否改善冷态起点。无 TTM、无电子系统、无外部热源；该设置仅用于预平衡筛查，不是生产边界条件。

## 设置

- 结构：917,378 原子 QSGS 模型，`boundary p p f`；起点沿用 Ag-only FIRE 200 步 Ag 坐标及原始骨架。
- 势函数：Ag EAM、Ti–Si–C TSCBOP/Tersoff 与 Morse calibration-v0，均未改参数。
- 底层 0–20 Å 固定；其余活动区使用 NVT 300 K，阻尼时间 0.1 ps；活动 Ti/Si/C 另加 `spring/self K=10 eV/Å²`。
- 时间步长 0.1 fs，运行 250 步（0.025 ps），4 线程 OpenMP。

## 结果

| 时间 (ps) | 活动区 (K) | Ag (K) | Ti/Si/C (K) | Ag 最大位移 (Å) | 骨架最大位移 (Å) |
|---:|---:|---:|---:|---:|---:|
| 0.000 | 300.00 | 299.43 | 300.17 | 0.000 | 0.000 |
| 0.005 | 326.89 | 300.35 | 334.62 | 0.040 | 0.121 |
| 0.010 | 355.68 | 303.70 | 370.82 | 0.078 | 0.220 |
| 0.015 | 328.02 | 310.72 | 333.05 | 0.116 | 0.333 |
| 0.020 | 298.56 | 322.95 | 291.45 | 0.156 | 0.442 |
| 0.025 | 339.94 | 340.91 | 339.65 | 0.202 | 0.538 |

运行完整结束，917,378 个原子均保留。移除 NVT 和弹簧后执行 `run 0`，活动 Ag 与骨架最大 pair-force 模分别约为 4.4438 和 4.2712 eV/Å。0.025 ps 仅为 0.1 ps 阻尼时间的四分之一，温度仍强烈振荡且末温约 340 K；不能把短时温度曲线解释为热平衡。

## 判断

全活动区 NVT 与弱骨架弹簧下，25 fs 内没有得到 300 K 平台，解除人工弹簧后的最大物理力仍较高。该试验说明单纯更换冷态预平衡控温分组没有解决当前初态的残余应力释放；它不证明长期一定发散，也不验证正式 TTM 条件。复合冷态基线仍未通过，暂不接入热源。

## 复现文件

- 输入：[`in.cold_pre_equil_K10_NVT_0p025ps`](in.cold_pre_equil_K10_NVT_0p025ps)。
- 日志：[`log.cold_pre_equil_K10_NVT_0p025ps_omp4`](log.cold_pre_equil_K10_NVT_0p025ps_omp4)，SHA-256 `5bea1a6817f27a6f964054d842f9e93778e1bce802afd81ccd4f96dd7c4dba1e`。
- 终点活动区坐标和速度：[`final_mobile_pre_equil_K10_NVT_0p025ps.dump.gz`](final_mobile_pre_equil_K10_NVT_0p025ps.dump.gz)，SHA-256 `e547c7cff77750ab7483de8aa0861d654b8efdc8a44263f93a38232371936d59`。

```bash
OMP_NUM_THREADS=4 lmp -sf omp \
  -in cold_pilot/in.cold_pre_equil_K10_NVT_0p025ps \
  -log cold_pilot/log.cold_pre_equil_K10_NVT_0p025ps_omp4
```
