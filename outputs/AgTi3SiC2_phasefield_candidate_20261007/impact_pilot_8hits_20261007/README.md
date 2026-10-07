# Ag–Ti₃SiC₂ 双连续小模型：8 发轰击熔坑候选

本目录是按实际入射点校正局部区域后重新运行的主结果。相场双连续掩膜、相组成和晶格映射保持不变。

- **结果报告：** [REPORT_中文.md](REPORT_中文.md)
- **形貌预览：** [surface_morphology_8hits_1p80ps.png](surface_morphology_8hits_1p80ps.png)
- **初始结构：** [structures/active_local_relaxed.data](structures/active_local_relaxed.data)
- **末态结构：** [structures/impact_8hits_final_cool_0p20ps.data](structures/impact_8hits_final_cool_0p20ps.data)
- **无轰击末态：** [structures/nohit_8shot_final_cool_0p20ps.data](structures/nohit_8shot_final_cool_0p20ps.data)
- **末段 OVITO 轨迹：** [dumps/impact_8hits_final_cool_0p20ps.dump](dumps/impact_8hits_final_cool_0p20ps.dump)
- **最终几何/温度审计：** [audit/eight_hit_final_audit.json](audit/eight_hit_final_audit.json)
- **逐段输入、日志、restart、dump 和结构：** 分别见 `inputs/`、`logs/`、`restart/`、`dumps/`、`structures/`。

8 次 62 eV Ag 入射，累计模拟 1.80 ps。末态中心 2 Å 内最高基底原子比同时间无轰击对照低约 11.0 Å，中心没有基底原子高于 35 Å，呈熔坑样凹陷；基底原子未丢失，两相主体连通。

这个结果是探索性熔坑候选，不是经验证的真实电弧熔池：末态仍处于高温非平衡状态，有一个 1.25 Å Ag–Ag 近接，且 Ag–Ti₃SiC₂ 交叉势未独立验证。见报告的限制说明。

`SHA256SUMS.txt` 列出本目录文件的 SHA-256，清单自身不包含在内。
