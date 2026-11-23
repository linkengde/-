# Ag/Ti3SiC2 固定尺寸同步加粗原子模型

## 当前版本：两相同步加粗，质量比保持 1:1

以仓库中的纳米原子模型为基准，在固定材料区和模拟盒内同步增加 Ti3SiC2 骨架与 Ag 网络。当前可查看版本是 `step05`。

- 材料区仍为 **200 × 200 × 200 Å（20 nm）**；模拟盒仍为 **200 × 200 × 300 Å（30 nm）**。
- 边界条件为 x/y 周期、z 非周期。
- Ti3SiC2 骨架界面向外扩 **0.80 Å/侧**；Ag 从原银网络附近的 Ag 晶格位点向外扩展。
- Ti3SiC2 与 Ag 的净原子数分别增加 **14,796** 和 **4,474**，约为原相原子数的 **4.41%**。
- 最终 Ti3SiC2 为 **350,514** 个原子，Ag 为 **105,991** 个原子，总计 **456,505** 个原子。
- Ag/Ti3SiC2 质量比为 **0.99999936**，按原子整数计数实现 1:1。
- 两相各自为单一连通网络；最小跨相原子距离为 **2.2000 Å**。

### 文件

- OVITO 可直接打开的完整原子模型：`output/step05_both_phases/ag_ti3sic2_both_phases_thickened_step05.extxyz`
- LAMMPS atomic 数据文件：`output/step05_both_phases/ag_ti3sic2_both_phases_thickened_step05.data`
- 相场和尺寸数据：`output/step05_both_phases/both_phases_step05_phase_field.npz`
- 建模与质量检查报告：`output/step05_both_phases/both_phases_step05_report.json`
- OVITO 原始模型/同步加粗模型对比图：`visualization/Ag_Ti3SiC2_both_phases_thickened_step05.png`
- 对比图出版格式：同名 `.pdf` 和 `.tif`
- 可重复建模和渲染程序：`build_step05_both_phases.py`、`render_step05_both_phases_ovito.py`

### 两相的扩展规则

骨架从原 QSGS 距离场向外拓宽 0.80 Å/侧。Ag 从原 Ag 网络附近的晶格位点扩展，并允许占据原 Ti 相几何掩膜中的原子级空隙；与 Ti、Si、C 原子的距离至少为 2.2 Å。随后从连通 Ag 主网络中取出达到 1:1 质量比的原子数，并重新检查连通性。

这个固定尺寸模型的加粗幅度受到可用空隙、Ag/Ti3SiC2 质量比和相间距离的共同限制。当前报告记录了 Ag 位点候选的连通分量、保留的基准 Ag 原子数和质量比修整数量，便于检查银相扩展是否足够连续。

## 历史版本

`step01`–`step03` 是未施加 1:1 质量比约束的早期加粗候选。`step04` 曾只从扩大骨架掩膜外选 Ag，因此虽然达到 1:1 质量比，却没有按当前要求让 Ag 网络也向空隙扩展。`step05` 是后续采用两相同步扩展规则的版本。

## 重建

从仓库根目录运行：

```bash
python models/ag_ti3sic2_qsgs_200A_skeleton_thickening/build_step05_both_phases.py
```

渲染脚本使用 OVITO 3.17.2 的 StandardRenderer：

```bash
python models/ag_ti3sic2_qsgs_200A_skeleton_thickening/render_step05_both_phases_ovito.py
```

## 限制

- 这是纳米尺度模型；骨架加粗不会将材料实际尺寸变成微米。
- 当前结构未做能量最小化、冷态稳定性检查或 TTM/MD 加热验证。
- 质量比、几何连通性和相间距离通过了脚本检查，不代表此前未通过的完整混合势力学稳定性已经解决。
