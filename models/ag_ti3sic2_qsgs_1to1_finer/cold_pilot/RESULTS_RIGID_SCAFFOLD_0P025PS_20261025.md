# 刚性骨架下 Ag 无源冷态筛查：0.025 ps

日期：2026-10-25（Asia/Shanghai）。本测试将 Ti₃SiC₂ 骨架完全固定，用来把骨架运动与 Ag 响应分开；它不是柔性骨架验证、平衡计算或熔池测试。

## 设置

- 917,378 原子 QSGS 模型，`boundary p p f`；使用与前期冷态筛查相同的 Ag-only FIRE 200 步坐标。
- 使用现有 Ag EAM、Ti–Si–C TSCBOP/Tersoff、Morse calibration-v0 势，参数未改。
- Ti/Si/C 与底层固定；活动 Ag 使用 NVT 300 K，温度阻尼时间 0.1 ps；无外部热源、无 TTM。
- 时间步长 0.1 fs，运行 250 步（0.025 ps），LAMMPS 22 Jul 2025 – Update 4，4 线程 OpenMP 势函数。

## 结果

| 时间 (ps) | Ag 温度 (K) | Ag MSD (Å²) | Ag 最大位移 (Å) | Ag 最大 pair-force 模 (eV/Å) |
|---:|---:|---:|---:|---:|
| 0.0000 | 300.00 | 0 | 0 | 3.8064 |
| 0.0050 | 301.03 | 1.7346×10⁻⁴ | 0.0414 | 3.8473 |
| 0.0100 | 304.15 | 6.9420×10⁻⁴ | 0.0821 | 3.8867 |
| 0.0150 | 309.79 | 1.5635×10⁻³ | 0.1220 | 3.9260 |
| 0.0200 | 318.64 | 2.7843×10⁻³ | 0.1610 | 3.9585 |
| 0.0250 | 331.59 | 4.3630×10⁻³ | 0.2038 | 3.9735 |

250 步完成，原子数保持 917,378，日志没有丢原子报告。Ag 温度随时间升高约 31.6 K，最大位移增至约 0.204 Å；末态 Ag 最大 pair-force 模约 3.9735 eV/Å。温度阻尼时间是模拟时长的 4 倍，因此该 25 fs 检查不足以证明 300 K 热平衡。虽然刚性骨架抑制了骨架运动，但 Ag 仍持续响应当前结构与界面势的残余力。

## 判断

完全固定骨架的 25 fs 短跑比 K=10 弱约束测试更少地把势能释放到骨架热运动中，但 Ag 仍升温，且该测试人为去掉了骨架的热响应。它可以保留为“刚性几何约束下的数值对照”，不能当作已验收的复合冷态起点，也不能证明真实 Ti₃SiC₂ 对熔融 Ag 的束缚或传热效果。当前复合体系冷态基线仍未通过；尚未运行 TTM 加热。

## 复现文件

- 输入：[`in.fixed_skeleton_ag_cold_0p025ps`](in.fixed_skeleton_ag_cold_0p025ps)。
- 日志：[`log.fixed_skeleton_ag_cold_0p025ps_omp4`](log.fixed_skeleton_ag_cold_0p025ps_omp4)，SHA-256 `78907cac3d16adc24d49637d20ad8437dc382d2ffabbee3ed83c85f1aa866d88`。
- 终点活动 Ag 坐标和速度：[`final_mobile_ag_rigid_scaffold_0p025ps.dump.gz`](final_mobile_ag_rigid_scaffold_0p025ps.dump.gz)，SHA-256 `117732cbf41e83099bbf2b1c469df730796b2499da63fa57060e77dc437a9c5c`。

```bash
OMP_NUM_THREADS=4 lmp -sf omp \
  -in cold_pilot/in.fixed_skeleton_ag_cold_0p025ps \
  -log cold_pilot/log.fixed_skeleton_ag_cold_0p025ps_omp4
```
