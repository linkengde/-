# 第466发 0.025 fs 时间步敏感性复核

本目录从第465发末态独立复制输入结构，重新执行第466发同一颗 62 eV Ag 入射。初态 data 的 SHA-256 与 100 ps 后轰击诊断分支的 `after_hit_465.data` 完全相同。位置、速度、ID、热浴参数、随机种子、势文件、边界和0.2 ps物理时长都保持一致；只把时间步从0.05 fs减半为0.025 fs，并把积分步数从4000加倍到8000。轨迹每2.5 fs输出一帧。

分析识别出 ID 4771 是本次入射 Ag，ID 4508 是基体 Ag。时间步减半后仍观察到相同碰撞：在两支路共同的103.905 ps采样点，原子间距都是约1.326515 Å，力约151.54 eV/Å。更密的采样在103.9025 ps记录到1.316834 Å、约163.34 eV/Å的瞬时近距和峰值力。说明 5 fs 原始采样漏掉了碰撞峰值；近接事件本身可重复，但其瞬时力峰不应当作收敛的定量指标。Ag–Ag pair 叠加 ZBL，是这个高能入射接触中的短程排斥响应。

两条积分的热浴能量校正总能量 `Etot + f_bathfix` 在整段0.2 ps内波动范围分别为0.00734 eV和0.00508 eV。没有丢原子、NaN/Inf或LAMMPS错误。末态 Ti₃SiC₂ 最大连通簇为2905/2908（99.90%），Ag最大连通簇为1793/1863（96.24%）。

- 复核报告：[`REPORT_dt_sensitivity_中文.md`](REPORT_dt_sensitivity_中文.md)
- 近接距离与力：[`images/dt_sensitivity_hit466.png`](images/dt_sensitivity_hit466.png)
- 可复算分析：[`analysis/analyze_dt_sensitivity.py`](analysis/analyze_dt_sensitivity.py)
- 逐帧配对数据：[`analysis/pair_force_timeseries.csv`](analysis/pair_force_timeseries.csv)
- 原始数据及输出：`inputs/`、`structures/`、`potentials/`、`logs/`、`dumps/`、`restart/`

这是一发事件的减半步长检查，不构成系统性时间步收敛研究，也不验证 Ag–Ti₃SiC₂ 交叉势。
