# 740 eV Ag 表层局部脉冲随机种子复算报告

日期：2026-11-25

## 目的

检验首个 739.62 eV 局部窄脉冲产生的浅凹和中心喷出是否对随机速度方向敏感。该复算从相同 1.700 ps 起始快照重新独立启动，只更换随机种子，并将速度扰动整体缩放到与首次计算相同的实际动能剂量。

## 完全匹配的设置

- Ag 数量：106,396；全体系：458,248 原子。
- 热源中心：(100,100,203) Å；横向高斯 σ=12 Å；向内指数深度=6 Å。
- 添加 Ag 动能：739.62464794 eV；首个种子与种子2剂量差小于数值舍入。
- 随机种子：20261130。
- LAMMPS `metal` 单位，NVE，时间步长0.00025 ps（0.25 fs），观察2.0 ps。
- Ag–Ag：Ag_u3.eam；骨架固定；探索性交叉作用和边界均与首个试算一致。没有恒温器、TTM或连续热源。

## 结果

- 末态仍有全部 458,248 原子；总能量全程范围约0.19 eV。
- 末态全 Ag 温度2371 K，热点温度3631 K，Ag 最大力14.24 eV/Å。
- 中心 r<20 Å 且 180≤z<220 Å 的 Ag 从573降至422；这是上表层稀疏/结构重排指标，不等于这些原子都已蒸发。
- 分箱表面高度中位数：初态中心/环带为201.96/200.84 Å；末态为209.77/208.50 Å。中心相对环带没有下陷。
- 热源 r<25 Å 内，轨迹中最多2个 Ag 同时高于z=250 Å；最高局部 z=265.32 Å。末态仍有2个处于该高度以上，但都仍在盒内。全场有8个 Ag 高于z=250 Å，其中包括远离热源的原子。

首个 739.62 eV 种子的中心表面相对环带末态低约2.27 Å；种子2末态反而高约1.27 Å。两次均有中心上层稀疏和少量喷出，但没有重现稳定中心凹陷。当前结果支持“局部喷射/稀疏化倾向”，不支持“已形成可重复熔坑”。

## 复现入口和文件

- `prepare_gaussian_pulse_surface_local_33750_seed2.py`：同剂量速度扰动，种子20261130。
- `prepared/local_gaussian_pulse_surface_local_33750_seed2_start.dump`：实际运行起始坐标和速度。
- `in.local_gaussian_pulse_surface_local_33750_seed2_2ps.lmp`：LAMMPS输入。
- `log.local_gaussian_pulse_surface_local_33750_seed2_2ps.lammps`：完整运行日志。
- `results/surface_local_33750_seed2_2ps_analysis.json`、`results/thermo_audit_surface_local_33750_seed2_2ps.json`、`results/trajectory_metrics_surface_local_33750_seed2_2ps.csv`：空间、能量和逐帧喷出统计。
- `trajectory_surface_local_33750_seed2_2ps/final_all_atoms_surface_local_33750_seed2.dump`：末态全原子快照。
- `results/center_slice_comparison_surface_local_33750_seed2_2ps_Chinese.png`：OVITO 中文截面对比图。

全体系 Ti₃SiC₂ 骨架在本路线中固定，因此该试算不能描述骨架吸热或晶格传热；局部脉冲为速度扰动而非 TTM 激光沉积。结果仅用于探索和筛选下一轮剂量，不作为已验证熔坑机制的证据。
