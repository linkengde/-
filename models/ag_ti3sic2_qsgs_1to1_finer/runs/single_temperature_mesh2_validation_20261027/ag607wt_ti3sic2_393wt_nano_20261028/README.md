# Ag 60.7 wt% / Ti₃SiC₂ 39.3 wt% 纳米几何预览

本目录只包含新配比几何、组成计算、形貌审计与预览图。**尚未启动新配比的 1× 或 2× 求解**，等几何审阅后再建模。新配比没有使用任何 Ag70 检查点。

## 配比换算

按给定称量和候选密度：

- Ag：60.7 g / 10.49 g·cm⁻³ = 5.7864633 cm³
- Ti₃SiC₂：39.3 g / 4.52 g·cm⁻³ = 8.6946903 cm³
- Ag 理论固相体积分数：0.3995857977（39.9585798 vol%）
- Ti₃SiC₂ 理论固相体积分数：0.6004142023（60.0414202 vol%）

换算采用 φᵢ=(mᵢ/ρᵢ)/Σⱼ(mⱼ/ρⱼ)。模型没有显式孔隙相，因此体积分数是无孔固相归一化结果；粉末压实孔隙、烧结收缩及两相实际密度未提供，无法从称量推算含孔总体积中的体积分数。

## 几何

几何由已存在的 QSGS 生长秩场阈值化，再沿用原交替最大面连通分量清理算法生成。随机生长没有重跑；原始 QSGS 父几何未修改。

- 原生网格：100×100×50
- 物理尺寸：299.635×301.857×158.6115 Å
- 原生格距：2.99635×3.01857×3.17223 Å
- Ag：199,793 voxels，39.9586 vol%
- Ti₃SiC₂：300,207 voxels，60.0414 vol%
- 离散体积分数误差：每相 2.0230×10⁻⁷，低于一个体素的分辨率
- 两相在 x/y 周期、z 开放边界下各为一个连通分量；x、y 均有周期绕行，且均贯通 z 的上下表面

两倍网格计划直接对这同一二值几何使用 solver 的 `mesh_refinement=2`，形成 200×200×100；不会重新随机生成形貌。

## 几何细节与限制

界面距离直径是基于约 3 Å 二值栅格的 voxel-center 欧氏距离描述，不是原子级最小颈宽。按周期 x/y 和开放 z 计算，Ag 中约 38.4% 的相内 voxel 距离直径小于 6 Å，约 78.4% 小于 10 Å；Ti₃SiC₂ 对应约 26.3% 和 62.1%。该几何有较细的栅格颈部，低于一个 voxel 的细节不可分辨，后续结果需保留这一限制。定向中位弦长见 `geometry_audit.json`。

## 复现

从仓库根目录依次运行：

```bash
python models/ag_ti3sic2_qsgs_1to1_finer/runs/single_temperature_mesh2_validation_20261027/ag607wt_ti3sic2_393wt_nano_20261028/build_geometry.py
python models/ag_ti3sic2_qsgs_1to1_finer/runs/single_temperature_mesh2_validation_20261027/ag607wt_ti3sic2_393wt_nano_20261028/audit_geometry.py
python models/ag_ti3sic2_qsgs_1to1_finer/runs/single_temperature_mesh2_validation_20261027/ag607wt_ti3sic2_393wt_nano_20261028/preview_geometry.py
```

`geometry_preview.png` 显示抽样三维相点云和 XY/XZ/YZ 三个中面；完整三维贯通由程序化面连通与绕行审计验证。检查摘要在 `composition_geometry.json`，完整形貌指标在 `geometry_audit.json`，文件哈希在 `SHA256SUMS.txt`。

## 输入来源与旧算例状态

- 原 QSGS 父几何 SHA256：`7e0e1509b414597ce4e2a1c9fc387ed9656caf507f1186e2d28de79c3c1bfc43 `
- 本次新几何 SHA256：`edebee1d8d254a60e51c7de0b63d96c10263020989c9c963449b932323e95695`
- 停止前 checkout 为 detached HEAD `940debd8eff1d293edc3543585c85b8bf4c74f40`，与目标分支当时的远端提交一致。Ag70 2× 任务已停止；其 1.31 ps 检查点恢复验证通过并已推送。旧目录、几何、结果和检查点均保留。

