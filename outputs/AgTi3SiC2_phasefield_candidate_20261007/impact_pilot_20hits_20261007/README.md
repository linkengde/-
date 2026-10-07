# Ag–Ti₃SiC₂ 双连续模型：20 发轰击后的熔坑样凹陷

本目录保存修正底部支撑边界后完成的 20 发 Ag 轰击小模型。坑形图已在对话中直接预览，也可在 `images/crater_20_final.png` 查看。

## 打开模型

- OVITO 直接加载 `structures/cool_05_quench_0p5ps.data` 查看最终原子结构。
- 加载 `dumps/impact_20_endframes_trajectory.dump` 可逐帧查看轰击及淬冷的阶段末态；每发只保留结束帧，时间标签写为累计飞秒，详见 `dumps/trajectory_index.csv`。
- 类型映射：type 1=Ti、2=Si、3=C、4=Ag。可将 Ti/Si/C 显示为紫色，将 Ag 显示为金色。
- `images/crater_20_final.png` 左图是由表面原子高度按半径统计后重建的径向平均坑形；右图是通过实际末态原子的竖直切片，紫色为 Ti₃SiC₂ 骨架、金色为 Ag。左图是形貌统计可视化，不是原子快照。

## 主要文件

- `structures/base_active_local_relaxed.data`：轰击前 4,305 原子模型。
- `structures/cool_03_0p5ps.data`：前 12 发后的阶段检查点，后续 13–20 发从此继续。
- `structures/cool_05_quench_0p5ps.data`：最终模型，4,325 原子。
- `restart/cool_05_quench_0p5ps.restart`：可供 LAMMPS 续算的最终 restart。
- `inputs/`、`logs/`：各段输入与日志；`potentials/` 保存本次使用的势文件。
- `audit/`：表面高度、相连通性、原子近邻及势文件哈希审计。
- `SHA256SUMS.txt`：本目录文件校验清单（清单自身不计入）。

完整数值与结论边界见 [中文结果报告](REPORT_中文.md)。
