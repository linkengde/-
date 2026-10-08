# 100 ps 末态后的液态判据诊断

从 `to100ps_200fs_20261008/checkpoints/after_hit_446` 复制状态，独立追加 25 个 Ag 入射（第 447–471 发），每发 0.2 ps，总计 5 ps。保留同一 62 eV 入射、200 fs 间隔、温控、边界条件和势文件；只加密轨迹输出，每 100 步（5 fs）写一帧，坐标同时记录 wrapped 与 unwrapped 形式，并记录速度和受力。

本目录的 `state_current` 和检查点属于诊断副本，不覆盖 100 ps 分支。每发完成后检查 LAMMPS 错误/丢原子、有限数值、最大力、原子数和 TSC/Ag 最大连通簇。任何严重事件会停在已保存状态。每发 dump 在审计后 gzip 压缩。

轨迹用于分析局部 Ag 的时间相关结构与扩散，并检查 TSC 骨架是否仍连通。温度、单帧形貌或单独的结构无序都不能证明熔化；Ag–TSC 交叉势尚未独立验证，结论仍是探索性模型结果。

## 诊断结果

轰击支路与匹配的无轰击支路都从同一 100 ps 状态开始，各保存 1001 帧、每 5 fs 一帧。轰击支路的 227 个起始坑区 Ag 示踪原子，漂移校正 MSD 在 1 ps 为 19.96 Å²、在 5 ps 为 99.55 Å²；对照分别为 5.44 和 31.14 Å²。两支路的局部平均 Q₆ 在 5 ps 内大致不变（轰击 0.3875→0.3883；对照 0.3875→0.3830）。所以数据支持“继续轰击提高了局部 Ag 迁移”，尚未形成独立、明确的液相判定。

5 fs 采样还发现第 466 发附近一次短暂 Ag–Ag 近接：原子 4508 与 4771 在 103.905 ps 相距 1.327 Å，单原子力约 151.5 eV/Å；前后 5 fs 距离约 1.500 Å、1.785 Å。该事件未导致丢原子或 LAMMPS 数值报错，但需要作时间步敏感性复核。

- 完整报告：[`REPORT_liquid_diagnostic_5ps_中文.md`](REPORT_liquid_diagnostic_5ps_中文.md)
- 轰击/对照图：[`images/impact_vs_control_5ps.png`](images/impact_vs_control_5ps.png)
- 匹配无轰击对照：[`../liquid_control_noimpact_5ps_20261008/README_zh.md`](../liquid_control_noimpact_5ps_20261008/README_zh.md)

主驱动：`run_liquid_diagnostic.py`。起始 data/restart 和势文件均从 100 ps 分支复制；无轰击支路使用同一初态。报告中的 MSD 拟合是非平衡“表观迁移率”，不应解释为平衡自扩散系数。
