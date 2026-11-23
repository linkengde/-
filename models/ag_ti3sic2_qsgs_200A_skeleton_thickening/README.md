# Ag/Ti3SiC2 固定尺寸同步加粗原子模型

## 当前版本：两相网络进一步粗化加粗，质量比保持 1:1

以仓库中的纳米原子模型为基准，在固定材料区和模拟盒内同步增加 Ti3SiC2 骨架与 Ag 网络。`step05` 通过占据空隙增加了两相原子数；`step06` 进一步平滑细枝、将相内原子重分配到较粗的连通枝干，同时保持 step05 的相原子数和质量比。

- 材料区仍为 **200 × 200 × 200 Å（20 nm）**；模拟盒仍为 **200 × 200 × 300 Å（30 nm）**。
- 边界条件为 x/y 周期、z 非周期。
- step05 将 Ti3SiC2 骨架界面向外扩 **0.80 Å/侧**，Ag 从原银网络附近的 Ag 晶格位点向外扩展；两相原子数相对原始纳米模型各增加约 **4.41%**。
- step06 使用 **1.5 个网格体素（约 3.75 Å）**的周期 x/y 高斯平滑重新分配相内原子；两相的总原子数保持不变，局部枝干更粗、细枝更平滑。
- 当前 step06：Ti3SiC2 **350,514** 个原子，Ag **105,991** 个原子，总计 **456,505** 个原子。
- Ag/Ti3SiC2 质量比为 **0.99999936**，按原子整数计数实现 1:1。
- 两相各自为单一连通网络；最小跨相原子距离为 **2.2000 Å**。
- 等体积网格掩膜的平均局部半厚度代理值：Ti3SiC2 从 **2.64 Å** 增至 **3.96 Å**；Ag 从 **2.53 Å** 增至 **3.46 Å**。这是几何厚度指标，不是经过弛豫的材料界面宽度。

### 文件

- 当前 step06 的 OVITO 原子模型：`output/step06_coarsened_networks/ag_ti3sic2_coarsened_networks_step06.extxyz`
- step06 LAMMPS 数据文件、相场、报告及校验和：`output/step06_coarsened_networks/`
- step05/step06 OVITO 对比图：`visualization/Ag_Ti3SiC2_network_coarsened_step06.png`，以及同名 `.pdf`、`.tif`
- step05 同步扩展版原子模型及报告：`output/step05_both_phases/`
- 原模型/step05 对比图：`visualization/Ag_Ti3SiC2_both_phases_thickened_step05.png`，以及同名 `.pdf`、`.tif`
- 可重复建模和渲染程序：`build_step05_both_phases.py`、`build_step06_coarsened_networks.py`、`render_step05_both_phases_ovito.py`、`render_step06_coarsened_networks_ovito.py`

### 两相的扩展规则

骨架从原 QSGS 距离场向外拓宽 0.80 Å/侧。Ag 从原 Ag 网络附近的晶格位点扩展，并允许占据原 Ti 相几何掩膜中的原子级空隙；与 Ti、Si、C 原子的距离至少为 2.2 Å。随后从连通 Ag 主网络中取出达到 1:1 质量比的原子数，并重新检查连通性。

这个固定尺寸模型的加粗幅度受到可用空隙、Ag/Ti3SiC2 质量比和相间距离的共同限制。当前报告记录了 Ag 位点候选的连通分量、保留的基准 Ag 原子数和质量比修整数量，便于检查银相扩展是否足够连续。

### step06 的继续加粗方式

固定盒子和 1:1 质量比时，继续增加两相原子数会很快耗尽可用 Ag 连通位点。step06 不改变相原子总数，而是平滑两相的空间占位密度，并在原有 Ti3SiC2 与 Ag 晶格上重新选位，使细枝减少、主枝变粗。网格半厚度统计保存在 `output/step06_coarsened_networks/coarsened_step06_report.json` 中。该步骤是未弛豫的几何重分配，需在用于 MD 前再做力学稳定性检查。

## 历史版本

`step01`–`step03` 是未施加 1:1 质量比约束的早期加粗候选。`step04` 曾只从扩大骨架掩膜外选 Ag，因此虽然达到 1:1 质量比，却没有按当前要求让 Ag 网络也向空隙扩展。`step05` 是后续采用两相同步扩展规则的版本。

## 重建

从仓库根目录运行：

```bash
python models/ag_ti3sic2_qsgs_200A_skeleton_thickening/build_step05_both_phases.py
python models/ag_ti3sic2_qsgs_200A_skeleton_thickening/build_step06_coarsened_networks.py
```

渲染脚本使用 OVITO 3.17.2 的 StandardRenderer：

```bash
python models/ag_ti3sic2_qsgs_200A_skeleton_thickening/render_step05_both_phases_ovito.py
python models/ag_ti3sic2_qsgs_200A_skeleton_thickening/render_step06_coarsened_networks_ovito.py
```

## 限制

- 这是纳米尺度模型；骨架加粗不会将材料实际尺寸变成微米。
- 当前结构未做能量最小化、冷态稳定性检查或 TTM/MD 加热验证。
- 质量比、几何连通性和相间距离通过了脚本检查，不代表此前未通过的完整混合势力学稳定性已经解决。
