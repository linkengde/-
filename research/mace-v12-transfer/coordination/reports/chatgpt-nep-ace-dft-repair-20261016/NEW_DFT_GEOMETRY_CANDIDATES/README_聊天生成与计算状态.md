# 新增10个DFT输入候选：本次聊天界面直接建模并上传GitHub

日期：2026-10-16。本文件和10个extxyz、审计CSV由**ChatGPT当前聊天界面**直接解析GitHub旧结构、实施刚体位移或静态小位移、计算周期近邻并提交独立分支。**没有任何新DFT自洽计算，也没有机器学习势训练、松弛、MD或TTM**。

## 结果

- 6个135原子Ag21/Ti57/Si19/C38完整界面：同一Si终止、同一结构父族、2个横向配准×最短跨相距离2.35/2.60/2.95 Å。使用2D三斜晶胞周期最小镜像求得目标最短距离，移动整层21个Ag，114个TSC原子完全保留。
- 2个Ag32 FCC固态微扰结构：坐标高斯静态位移每轴标准差0.04/0.08 Å，按全部原子去质心偏移，种子固定；是**静态畸变候选，不是液体结构、不是300K/高温轨迹**。
- 2个Ti3SiC2体相48原子微扰结构：每轴标准差0.025/0.05 Å；仍是**静态畸变候选，并非热平衡高温MAX相**。
- 所有输出extxyz只包含元素、坐标、父原子序号和晶胞/PBC；无E、F、stress标签。每项最短近邻/接触信息在GEOMETRY_NUMERICAL_AUDIT.csv。

## 结构源与证据

原始仓库 `linkengde/-`，参考commit `4ee3daab860a846fec4ef0c0ca4b031a36cfc87e`：

1. 135原子源：`research/mace-v12-transfer/coordination/reports/window-e/AG_TISC_C_DFT_HANDOFF_20261015/existing_dft/135atom_full_relaxation/optimized_structure.extxyz`，SHA-256 `eda2586c63c6c3c793ffce99417e7c04569941b0a2a916fda1b5402b7ba52d9f`。
2. 32原子Ag源：`research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/inputs/Ag_baseline_PROPOSAL.extxyz`，历史几何源 SHA-256 `21138a56bc3138106c102c4eb207ea76e2060bb5315f454a937f3ce26308318f`。
3. 48原子Ti3SiC2源：`research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/inputs/Ti3SiC2_baseline_PROPOSAL.extxyz`，历史几何源 SHA-256 `157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4`。

## 物理与数值限制

- 135原子界面六点共用单一结构父族，只能作为同族训练/探索候选，**绝不能作为六个独立测试家族**；缺其他终止面、独立配准父族和液态Ag–MAX接触标签。
- Ag21分相是奇电子中性体系；DFT分相扣能须保留一致的自旋/占据和电子设置，不得与native total/free energy混用。
- 当前本聊天计算容器无GPAW、ASE、LAMMPS，约4.7 GiB可用内存，无MPI安装；**无法在本容器完成GPAW单点或验证39,486原子MD性能**。
- 真正DFT作业需另有经过核验的GPAW26.7/ASE3.29/PAW setup环境，先进行k-point、展宽、真空/dipole和力误差检查；不把原Γ界面结果当作完成数值收敛。
- 原始仓库数据未改动，产生文件来源标记为CHATGPT_CURRENT_CONVERSATION_DERIVED，**非E窗口新完成DFT**。

## 当前门槛
```
GEOMETRY_GENERATION_COMPLETED=YES
STATIC_NEIGHBOR_AUDIT=PASS_CANDIDATE_ONLY
NEW_DFT_LABELS_GENERATED=NO
LIQUID_AG_DFT_DATA_READY=NO
HIGH_T_TSC_DFT_DATA_READY=NO
NEP_ACE_TRAINING_ALLOWED=NO
```
