# Ag–Ti₃SiC₂ 更细枝状骨架 QSGS 几何候选（seed 20261141）

这是根据当前要求重新生成的模型：先用随机四参数 QSGS 生成 Ti₃SiC₂ 骨架，再以 Ag 晶格位点填充互穿空间。它是全新随机种子，不沿用前两版的骨架掩膜。

## 模型尺寸与分层

- 材料：200 × 200 × 200 Å
- LAMMPS 计算盒：200 × 200 × 500 Å；z=200–500 Å 为真空观察区
- 原子数：469,040
- 元素计数：Ti 179,214；Si 59,738；C 119,476；Ag 110,612
- Ag/Ti₃SiC₂ 质量比：1.02056:1（接近 1:1）
- z=0–20 Å：固定层候选，仅标记骨架原子
- z=20–40 Å：恒温层候选，仅标记骨架原子
- z=40–200 Å：牛顿层
- z=200–500 Å：真空区

模型文件中的 `pbc=[True, True, False]`，即 x/y 周期、z 非周期。分层数组是供后续 LAMMPS 输入使用的标记；尚未运行动力学。

## QSGS 设置与几何筛查

- 网格间距：1 Å；随机种子：20261141
- 四参数概率：成核 0.0016、前向生长 0.070、分枝 0.030、反向生长 0.00060
- QSGS 生长目标体素占比 0.84；最终骨架掩膜阈值 0.52
- 骨架分枝最大直径估算：7.48 Å
- Ag 与骨架之间采用 2.0 Å 几何筛选间距；Ag–Ag 连通筛选截断 3.15 Å
- Ti₃SiC₂ 与 Ag 的主网络均在几何图论筛查中贯通 x、y、z；每个 20 Å 立方分区均有原子

以上是几何构建及连通性筛查结果，不是势函数稳定性或物理性质验证。模型尚未弛豫，尚未进行 MD；开始模拟前仍需用选定势函数检查界面受力、近邻距离和结构稳定性。OVITO 图片中的球半径仅用于显示，未改变原子坐标或结构文件。

## 文件

- `qsgs_thin_1A_seed41_200x200x500.data`：LAMMPS atomic data，可直接导入 OVITO 检查
- `qsgs_thin_1A_seed41_200x200x500.extxyz`：扩展 XYZ，可直接导入 OVITO
- `AgTi3SiC2_更细骨架_seed41_OVITO_CN_300dpi.png`：中文全模型预览图
- `qsgs_thin_1A_seed41_mask.npz`、`qsgs_rank_1A_seed20261141.npy`：掩膜和 QSGS 排名场
- `qsgs_thin_1A_seed41_report.json`、`qsgs_1A_screen_seed41_results.json`：模型参数、计数和几何筛查记录
- `build_thinner_1A_seed41_model.py`、`generate_and_screen_seed41.py`、`render_seed41_ovito.py`：生成、筛查和渲染脚本
- `support/build_microtwig_candidate.py`：构建脚本使用的通用几何辅助函数
- `SHA256SUMS.txt`：本目录文件校验值（不含校验文件自身）

## 重建说明

在仓库根目录运行脚本。Python 环境需安装 NumPy、SciPy、ASE、Pillow 和 OVITO Python 模块。完整重建原子结构时，辅助函数还会读取仓库中的 Ti₃SiC₂ CIF 和既有结构参考数据：

- `research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif`
- `models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/relaxed_backbone_candidate/`（辅助几何数据）

由已保存的排名场重建 data/extxyz：

```bash
python models/target_200x200x300_qsgs_material200/v14_exploratory_hybrid_ag_20261124/qsgs_thinner_1A_seed41_rebuild_20261126/build_thinner_1A_seed41_model.py
```

重新生成随机排名场并执行候选参数筛查（耗时和内存占用较高）：

```bash
python models/target_200x200x300_qsgs_material200/v14_exploratory_hybrid_ag_20261124/qsgs_thinner_1A_seed41_rebuild_20261126/generate_and_screen_seed41.py
```

OVITO Python 模块可用时，生成中文预览图：

```bash
python models/target_200x200x300_qsgs_material200/v14_exploratory_hybrid_ag_20261124/qsgs_thinner_1A_seed41_rebuild_20261126/render_seed41_ovito.py
```

## 加热结果

seed41 Ag 均匀 NVT 加热记录、输入、日志、Q6/MSD 轨迹、OVITO PTM 逐帧统计和末帧见 `runs/seed41_heating_q6_ptm_msd_20261126/HEATING_REPORT_CN.md`。最终时间为 1.800 ps。加热采用 Ag EAM + 固定骨架几何障碍的探索性模型，不是 TTM。
