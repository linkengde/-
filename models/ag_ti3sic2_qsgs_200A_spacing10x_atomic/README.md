# Ag/Ti3SiC2 原子模型：固定盒尺寸、枝间尺度放大 10 倍

这个版本按“先生成完整双相原子模型，再从最终模型提取 Ti3SiC2 骨架和局部枝段”的顺序生成。它不是把先前的体素预览当作原子模型。

## 建模方法

1. 核对已有纳米基准模型：437,235 个原子；材料区 200 × 200 × 200 Å；模拟盒 200 × 200 × 300 Å；x/y 周期、z 非周期。
2. 从该纳米模型的 QSGS Ti3SiC2 相标签场截取 8 × 8 × 8 个体素，即 20 × 20 × 20 Å 局部形貌；将该形貌映射到仍为 200 × 200 × 200 Å 的材料区，空间倍率为 10。
3. 按原模型的晶格参数和取向重新放置 Ti3SiC2 晶格与 Ag FCC 晶格，并保留原模型的 Ag:Ti3SiC2 等质量比。材料区、模拟盒和各层 z 范围不变。
4. 从新生成的完整原子模型读取 Ti/Si/C 原子坐标，重建相场并计算骨架中心线，最后从同一骨架提取局部区域。

“10 倍”只表示固定 20 nm 材料区内的枝间形貌尺度放大；它没有把模型的物理盒子扩展到微米。这是保持原子数可计算、用于比较枝间距变化的粗化原子模型。

## 本次结构数据

- 原子总数：407,009。
- Ti3SiC2：312,510 个（Ti 156,255；Si 52,085；C 104,170）。
- Ag：94,499 个；Ag/Ti3SiC2 质量比约 1.000。
- Ti3SiC2 与 Ag 原子网络各为一个连通分量。
- 从原子位置重建后提取的中心线：1,097 个体素；显示管半径 7.5 Å，仅供读图，不是原子尺寸。
- 从最终 Ti 原子位置重建相表面后，局部两枝表面间距为 93.702 Å（9.370 nm），约为原模型参考间距 0.93719 nm 的 10.00 倍。
- 材料区：200 × 200 × 200 Å（20 nm）；模拟盒：200 × 200 × 300 Å（20 × 20 × 30 nm）。
- 固定层 0–20 Å；恒温层 20–40 Å；牛顿层 40–200 Å；溅射区 200–250 Å；蒸发区 250–300 Å。

模型沿用基准的等质量配比，因此 Ag 没有填满全部互补几何体素；按各晶格位点密度估算，Ti3SiC2、Ag 和未占据体积分别约占材料盒的 45.8%、20.1% 和 34.1%。这是基准配比对应的孔隙，不应解释为无孔致密材料。

## 文件

### 原子模型

- `output/ag_ti3sic2_spacing10x.extxyz`：推荐用 OVITO 打开；保留元素、相标签和各 z 层标签。
- `output/ag_ti3sic2_spacing10x.data`：LAMMPS atomic 格式坐标与质量。
- `output/spacing10x_atomic_phase_field.npz`：生成完整原子模型时使用的体素相场。
- `output/atomic_model_report.json`：模型尺寸、原子数、配比、连接性和来源记录。

### 从完整原子模型提取的骨架

- `output/Ti3SiC2_skeleton_from_atoms.stl`：从原子位置重建的完整 Ti3SiC2 骨架相表面，坐标单位 Å。
- `output/Ti3SiC2_local_branch_from_atoms.stl`：同一原子重建相表面中截取的 13 nm 局部枝段，坐标单位 Å。
- `output/local_branch_atoms.extxyz`：同一局部区域的原子裁剪，含 Ti3SiC2 和 Ag。
- `output/atomic_skeleton_extraction.npz`、`output/skeleton_extraction_report.json`：相场、中轴、局部框和实测枝间距。
- `visualization/Ti3SiC2_spacing10x_atomic_whole_and_local.png`、`.pdf`：整体骨架、局部枝段及 9.370 nm 实际枝间距标注。

## OVITO

在 OVITO 中选择 **File → Load File**，打开 `output/ag_ti3sic2_spacing10x.extxyz` 查看原子模型。打开 STL 文件可单独查看提取的整体骨架或局部枝段。STL 是骨架表面模型，不包含元素原子。

## 适用范围

该结构尚未做能量最小化、冷态稳定性检查或 TTM/MD 加热。现有 Ag–Ti/Si/C 界面势也未因此获得验证。它用于查看十倍枝间距粗化后的原子级几何，不应直接作为已验证的力学模型或微米尺度实物。

可从仓库根目录依次运行：

```bash
python models/ag_ti3sic2_qsgs_200A_spacing10x_atomic/build_atomic_model.py
python models/ag_ti3sic2_qsgs_200A_spacing10x_atomic/extract_skeleton_from_atoms.py
```
