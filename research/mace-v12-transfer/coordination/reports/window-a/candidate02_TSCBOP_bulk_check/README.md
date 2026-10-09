# 候选02 TSCBOP 固定骨架 DFT 对照

在已收敛 Ti₃SiC₂ 同一固定48原子结构上，用原始TSCBOP.tersoff（LAMMPS pair_coeff type order Ti Si C）做单线程run0，与k4/k6 GPAW有限展宽力比较。无弛豫、升温、运动或参数改动。按LAMMPS ID还原原子顺序，确认输入坐标一致。势和DFT输出哈希在report.json。

TSCBOP对两种网格均给出相同约−311.699856eV势能，但不比较未经零点对齐的绝对能量。全原子力向量RMSE约1.99eV/Å，最大单原子误差约2.71eV/Å；C约2.13、Ti约2.21、Si约0.0033eV/Å。误差稳定重现于两网格。对称体相小力无法解释Ti/C较大的力失配，说明候选02骨架项在此未弛豫 CIF 几何上本身就与DFT力不相容；不能指望只调Ag交叉项修复所有框架原子分量。

这仍是单一父结构的静态诊断，不说明热稳定性或已确认Stage69配置，也未拟合/采用新势。应先核验坐标/结构来源及该体相是否需要弛豫参考，再用体积、位移、剪切控制检查TSCBOP力响应；若失配持续，统一多元素拟合需重新覆盖Ti–C/Si骨架项。

复现：独立LAMMPS环境设LD_LIBRARY_PATH=/workspace/.venvs/lammps-cross/lib、OMP_NUM_THREADS=1、OPENBLAS_NUM_THREADS=1运行check.py。副本重跑以保留发布快照。
