# Ag/Ti3SiC2 原子模型：骨架加粗迭代

本目录从现有 437,235 原子纳米模型出发，对 Ti3SiC2 骨架做第一档小幅加粗，供查看和决定后续是否继续加粗。原模型没有覆盖或改写。

一并保留了建模所用的原始原子文件与相场输入：

- `models/ag_ti3sic2_qsgs_200A_ttm_finer/ag_ti3sic2_qsgs.extxyz`
- `models/ag_ti3sic2_qsgs_200A_ttm_finer/qsgs_morphology.npz`
- `models/ag_ti3sic2_qsgs_200A_ttm_finer/geometry_report.json`
- `models/ag_ti3sic2_qsgs_200A_ttm_finer/build_model.py`

## 按 1:1 质量比约束生成的版本

第一至第三档是早先用于观察加粗趋势的几何候选，Ag/Ti3SiC2 质量比发生了变化；它们不符合后来提出的 1:1 约束。`step04` 是按该约束重新生成的候选，建议按此版本继续比较。

- 在材料区和盒子尺寸不变、沿用当前两相晶格密度与界面间距时，最大测试可行加粗量为 **0.465625 Å/侧**。到 0.46875 Å/侧时，连通 Ag 位点数已少于 1:1 所需数量。
- `step04` 总原子数 **447,144**：Ti3SiC2 **343,326**，Ag **103,818**；计算 Ag/Ti3SiC2 质量比 **1.000005**，即按原子整数计数实现 1:1。
- 材料区为 **200 × 200 × 200 Å**，模拟盒为 **200 × 200 × 300 Å**，x/y 周期、z 非周期。
- 对比图：`visualization/Ag_Ti3SiC2_skeleton_thickened_equal_mass_1to1.png`。
- OVITO 模型：`output/step04/ag_ti3sic2_skeleton_thickened_step04.extxyz`。
- 结构报告：`output/step04/thickened_step04_report.json`。

在当前盒子和近似体相晶格密度设定下，继续明显加粗会挤掉达到 1:1 所需的 Ag 连通位点。因此 1:1 约束下的骨架增粗幅度很小；此前 1–3 档只是质量比变化的几何示例。

## 本版变化

- Ti3SiC2 相的界面向 Ag 一侧外扩 **1.0 Å/侧**，用原有晶格位点生成新增骨架原子。
- Ag 仅占据扩张后余下的可用位点，并保留一个连通 Ag 主相。
- 材料区仍为 **200 × 200 × 200 Å（20 nm）**；模拟盒仍为 **200 × 200 × 300 Å（30 nm）**；x/y 周期、z 非周期。
- 基准模型：Ti3SiC2 335,718 个、Ag 101,517 个，共 437,235 个原子。
- 加粗版：Ti3SiC2 356,166 个、Ag 101,472 个，共 457,638 个原子；Ag/Ti3SiC2 质量比约 0.942。

这只是第一档几何加粗。盒子物理长度仍为 20 nm；是否把后续达到的粗化形貌作为“微米风格模型”，由你看图后决定。它不会因为骨架变粗就自动变成真实微米尺寸。

## 第二档：累计外扩 2 Å/侧

第二档在第一档基础上继续加粗 1 Å/侧，仍从同一原始几何重建。第一档的 356,166 个 Ti3SiC2 原子全部保留；第二档新增 56,232 个 Ti3SiC2 原子，并移除余下区域中 18,433 个 Ag 原子。Ag 主相仍保持连通。

- 第二档总原子数：495,437；Ti3SiC2：412,398；Ag：83,039。
- Ag/Ti3SiC2 质量比约 0.666，低于第一档的 0.942。
- 材料区和模拟盒的物理尺寸、边界条件不变。
- 对比图：`visualization/Ag_Ti3SiC2_skeleton_thickened_step01_vs_step02.png`。
- 第二档 OVITO 原子文件：`output/step02/ag_ti3sic2_skeleton_thickened_step02.extxyz`。
- 第二档结构报告：`output/step02/thickened_step02_report.json`。

## 第三档：累计外扩 3 Å/侧

第三档再比第二档加粗 1 Å/侧，Ag 相仍保持连通并跨越 x/y 周期边界、贯穿 z 方向。

- 第三档总原子数：537,561；Ti3SiC2：490,554；Ag：47,007。
- 相对第二档增加 78,156 个 Ti3SiC2 原子，减少 36,032 个 Ag 原子。
- Ag/Ti3SiC2 质量比约 0.317，银相占比已明显下降。
- 对比图：`visualization/Ag_Ti3SiC2_skeleton_thickened_step02_vs_step03.png`。
- 第三档 OVITO 原子文件：`output/step03/ag_ti3sic2_skeleton_thickened_step03.extxyz`。
- 第三档结构报告：`output/step03/thickened_step03_report.json`。

## 查看

- 优先查看 `visualization/Ag_Ti3SiC2_skeleton_thickened_equal_mass_1to1.png`：原始模型与满足 1:1 质量比的加粗版 OVITO 对比图。
- 1:1 质量比版本的模型和报告位于 `output/step04/`。
- `visualization/Ag_Ti3SiC2_skeleton_thickened_step01_comparison.png`：原模型与加粗版的 OVITO 标准渲染并排图。
- `visualization/Ag_Ti3SiC2_skeleton_thickened_step01_comparison.pdf`、`.tif`：排版和打印版本。
- `output/ag_ti3sic2_skeleton_thickened_step01.extxyz`：OVITO 可直接打开的完整原子模型，含相和 z 层属性。
- `output/ag_ti3sic2_skeleton_thickened_step01.data`：LAMMPS atomic 数据文件。
- `output/thickened_step01_report.json`：尺寸、计数、连通性清理和限制。

## 重建

从仓库根目录运行（需要 Python、NumPy、SciPy 和 ASE）：

```bash
python models/ag_ti3sic2_qsgs_200A_skeleton_thickening/build_step01.py
python models/ag_ti3sic2_qsgs_200A_skeleton_thickening/build_step02.py
python models/ag_ti3sic2_qsgs_200A_skeleton_thickening/build_step03.py
python models/ag_ti3sic2_qsgs_200A_skeleton_thickening/build_step04_equal_mass.py
```

图由 OVITO 3.17.2 的 StandardRenderer 输出。当前云机无系统级 Vulkan 驱动，因此使用下载到 `/tmp` 的 Mesa lavapipe CPU Vulkan ICD 生成；该临时驱动不属于模型文件。

## 限制

- 模型尚未做能量最小化、冷态稳定性检查或 TTM/MD 加热。
- 前三档的 Ag/Ti3SiC2 质量比分别约 0.942、0.666、0.317，不满足 1:1；step04 按原子整数数目恢复到 1.000005。
- 本版用于看形貌和继续迭代，不代表势函数或界面力学已通过验证。
