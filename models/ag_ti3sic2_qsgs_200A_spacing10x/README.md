# 固定盒尺寸的枝间尺度粗化骨架模型

> **版本说明：本目录是体素相场几何预览，不含 Ag、Ti、Si、C 原子模型。** 之前的图和网格文件仅用于检查形貌；用户确认的流程要求完整原子模型先建立、随后从原子模型提取骨架。请查看 [原子级版本](../ag_ti3sic2_qsgs_200A_spacing10x_atomic/README.md) 和 `ag_ti3sic2_qsgs_200A_spacing10x_atomic/output/ag_ti3sic2_spacing10x.extxyz`。

本目录保留体素几何版本，便于追溯此前的形貌预览；它不能替代上面的原子级版本。

处理顺序是：先由原 QSGS 双相场取样并构建完整双相模型，再从完整 Ti₃SiC₂ 相提取中心线骨架，最后从同一副骨架截取局部分叉。材料区保持 **200 × 200 × 200 Å（20 nm）**；模拟盒保持 **200 × 200 × 300 Å（20 × 20 × 30 nm）**，名义边界条件为 x/y 周期、z 非周期。

尺度粗化由原 QSGS 完整场中的一个 20 Å 立方窗映射到 200 Å 材料区完成，空间倍率为 10。抽取后的骨架中心线用半径 7.5 Å 的管状表面显示。图中局部枝间距在该显示骨架的平滑表面上实测约 **9.414 nm**；原模型对应间距为 **0.93719 nm**。

## 文件

- `output/spacing10x_phase_field.npz`：完整两相标签、提取的中心线骨架和骨架显示管，网格 80 × 80 × 120，体素间距 2.5 Å。
- `output/spacing10x_composite.vtk`、`spacing10x_model_ti3sic2.vtk`、`spacing10x_model_ag.vtk`、`spacing10x_extracted_skeleton.vtk`：ASCII VTK 体素场，供 ParaView/VTK 使用；坐标单位为 Å。
- `output/spacing10x_model_ti3sic2.stl`、`spacing10x_model_ag.stl`：完整模型的两相独立表面；在 OVITO 中分别导入可单独着色。坐标单位为 Å。
- `output/spacing10x_ti3sic2_skeleton.stl`：完整提取骨架的表面网格。
- `output/spacing10x_local_skeleton.stl`：从完整骨架体素场中截取的 13 nm 局部分叉，可单独在 OVITO 中查看。
- `output/topology_metrics.json`：模型尺寸、体积分数、连通性、中心线节点和枝间距测量。
- `visualization/Ti3SiC2_spacing10x_whole_and_local.png`、`.pdf`：整体骨架及从同一模型截取的局部分叉图。
- `build_model_then_skeleton.py`、`spacing10x_geometry_utils.py`：模型生成和骨架提取程序；原 QSGS 输入为 `../ag_ti3sic2_qsgs_200A_ttm_finer/qsgs_morphology.npz`。

从仓库根目录安装依赖并复现：

```bash
python -m pip install -r models/ag_ti3sic2_qsgs_200A_spacing10x/requirements.txt
python models/ag_ti3sic2_qsgs_200A_spacing10x/build_model_then_skeleton.py
```

## OVITO

在 OVITO 中选择 **File → Load File**，导入 `spacing10x_model_ti3sic2.stl` 与 `spacing10x_model_ag.stl` 查看完整双相几何，导入 `spacing10x_ti3sic2_skeleton.stl` 查看整体骨架，或导入 `spacing10x_local_skeleton.stl` 查看局部分叉。

## 适用范围

这里的“10 倍”指枝间形貌尺度粗化。材料区和模拟盒仍是纳米尺寸。输出用于几何形貌展示，不含 Ag、Ti、Si、C 原子坐标；骨架显示管径也不是原子级结构尺寸，不能直接作为 LAMMPS 原子模型运行。x/y 周期是几何模型的名义边界设置，尚未验证原子级周期接缝或力学稳定性。
