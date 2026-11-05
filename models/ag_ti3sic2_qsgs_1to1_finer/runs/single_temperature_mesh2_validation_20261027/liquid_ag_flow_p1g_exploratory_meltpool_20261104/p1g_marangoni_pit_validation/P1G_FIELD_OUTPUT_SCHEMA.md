# P1G 三维场输出说明

## 原始快照和 HDF5

求解器每帧写一个原子替换的紧凑二进制快照，默认间隔 0.05 ps。完整重启检查点每 0.2 ps 保存一次，并在关热、完全凝固时额外保存；由检查点恢复后继续沿下一个 0.2 ps 边界保存。后处理器把固定网格快照转换成压缩 HDF5，字段维序为 time,x,y,z，网格为 64×64×64，记录每帧时间戳、等距网格坐标、边界周期方向和 Marangoni 开关。检查点与场快照分别保存。

求解器快照字段：

| 字段 | 单位 | 来源 |
|---|---|---|
| VOF vof | 1 | CLSVOF 样品/气体体积分数 |
| level set level_set_m | m | CLSVOF 距离场 |
| ag_marker, skeleton_marker | 1 | Ag 和 Ti₃SiC₂ 掩膜示踪量 |
| liquid_fraction | 1 | Ag 液相率 |
| temperature_K | K | 单温场 |
| enthalpy_density_J_m3 | J/m³ | 混合体积焓 |
| velocity_x/y/z_m_s | m/s | 速度场 |
| pressure_Pa | Pa | 压力场 |
| curvature_m-1 | 1/m | VOF 曲率 |
| surface_tension_candidate_N_m | N/m | 候选 Ag 表面张力 |
| vof_weighted_surface_tension_N_m | N/m | VOF 加权表面张力 |
| marangoni_acceleration_candidate_x/y/z_m_s2 | m/s² | 候选切向加速度 |
| source_power_density_W_m3 | W/m³ | 热源密度 |
| skeleton_and_mushy_drag_s-1 | 1/s | 骨架与糊状区阻力 |
| material_id | 1 | 0 气体、1 Ag、2 Ti₃SiC₂ |

HDF5 后处理另外存储由 level set 差分得到的界面法向量、由速度差分得到的涡量、每个 x-y 柱中最高 `VOF > 0.5` 网格单元对应的上表面包络高度，以及该位置的表面温度。没有样品的列记为 NaN。凹坑/凸起审计相对于初始包络高度计算，并要求位移超过半个网格间距；该高度仅适用于单值上表面，不能表示内部孔道气液界面或悬垂自由表面。后处理器已改为显式单位字典，避免把焓、表面张力和加速度写成错误单位。

## 派生表面量的解释

`free_surface_height_m` 使用每个 x-y 柱中最高的 `VOF > 0.5` 单元得到上表面包络高度；它只在顶面可表示为单值高度时适合解释为表面形状。凹坑深度、面积和体积只累计超过 `0.5 dx` 的列位移。法向量和涡量是后处理差分估算值，不是独立界面求解器直接输出的观测量。表面张力和 Marangoni 加速度标记为候选量，因为材料参数存在实验范围与外推限制。

## 同步和检查

快照头保存统一物理时间、网格间距和坐标原点；HDF5 数据集按时间帧存储。后处理会拒绝缺字段、截断快照、NaN/Inf、网格改变或时间戳/维度不一致的输入。A/B 0.01 ps 短测均输出 t=0 与 t=0.01 ps 两帧，文件大小各约 2.61 MB。
