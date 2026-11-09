# Round7 五帧完整体系轨迹与临时熔化候选筛查

## 可播放的完整体系轨迹

`round7_5frames_full_system_0p4_to_1p4ps.dump.gz` 是一个包含五个完整 LAMMPS dump 帧的压缩轨迹。帧时间为 **0.40、0.65、0.90、1.15、1.40 ps**；每帧均含全部 **352,637 个原子**。已用 OVITO 实际打开核对，识别为五帧。

下载该文件并用 OVITO 打开，在动画控制区点击播放即可查看五个时间点之间的变化。最后一帧为 1.40 ps；本次没有包含 1.50 ps。

## 五帧 Q6/PTM/MSD 临时筛查

同一组 6,448 个初始热源核心区域 Ag 原子在五帧上进行筛查。单原子候选条件为：窗口平均 Q6 < 0.408106928、五帧中至少 80% 为 PTM 非 FCC、1 ps 去整体漂移净位移平方 MSD ≥ 0.517167930 Å²。团簇以 3.5 Å Ag–Ag 连通距离识别，并考虑 X/Y 周期边界。

- 联合候选 Ag：2,667 个，占该热源核心区初始 Ag 的约 41.4%。
- 团簇：159 个；最大候选团簇 1,454 个 Ag。
- 最大团簇周期性包围范围约 83.5 × 171.3 × 33.6 Å。

这些阈值借用了 Round5 的纯 Ag 开发势探索性参考，尚未用 Round7 势重新标定。Round7 势也未独立验证液态 Ag 与界面行为，所以此结果是**熔化候选筛查**，不是已验证的熔化判定或材料性能结论。1.40 ps 检查点最大原子力约 159 eV/Å，且有一个原子进入 z>200 Å 区域；应谨慎解释。

## 文件

- `round7_5frames_full_system_0p4_to_1p4ps.dump.gz`：五帧完整体系、可在 OVITO 播放的压缩 LAMMPS 轨迹。
- `round7_1p4ps_joint_candidate_3d.png`：完整体系中候选 Ag 的三维示意图。
- `round7_1p4ps_joint_screen.json`：筛查规则、统计、团簇及限制。
- `round7_1p4ps_source_core_ag.csv`：热源核心区每个 Ag 原子的 Q6、PTM、MSD 和候选标记。
- `round7_1p4ps_joint_candidates.dump`：所有联合筛查候选 Ag。
- `round7_1p4ps_largest_candidate_cluster.dump`：最大候选团簇。
