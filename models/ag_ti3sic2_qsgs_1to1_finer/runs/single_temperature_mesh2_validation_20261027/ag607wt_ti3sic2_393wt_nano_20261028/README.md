# Ag 60.7 wt% / Ti₃SiC₂ 39.3 wt% 纳米几何预览

本目录包含新配比几何、组成计算、形貌审计、数值验收及 1×/2× 单温焓法算例。1 fs 短测已通过并逐位复现；完整计算从本几何和 300 K 初态开始，不使用 Ag70 检查点。

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

两倍网格直接对这同一二值几何使用 solver 的 `mesh_refinement=2`，形成 200×200×100；没有重新随机生成形貌。

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

## 热计算验收与完整算例

1× 与 2× 分别通过了 0.001 ps（1 fs）Case A 短测：源功率均为 `8.47205563778101×10⁻⁴ W`，能量相对误差绝对值低于 `7×10⁻¹⁶`，两次独立短测逐位一致，完整状态检查点恢复验证通过。短测结果在 `preflight/1x/preflight.json`、`preflight/2x/preflight.json` 和 `preflight_summary.json`。

完整算例复用现有求解器与物性、边界和热源参数。默认不设加热时间上限，最高温度保护为 2300 K；加热持续至稳定熔池判据达成并关闭热源，或触及温度保护。冷却阶段最多 5 ps。每 0.05 ps 保存完整检查点，液态连通与团簇量每 0.01 ps 采样。稳定熔池判据为：同一最大液态连通团簇连续存在至少 3 ps，且最新 2 ps 的团簇体积极差不超过均值的 5%。若触及温度保护，报告实际停机原因，不伪称形成稳定熔池。

从本目录启动两种网格的短测门控和并行完整算例：

```bash
python run_ag607wt_pair.py --max-cooling-ps 5 --max-temperature-K 2300
```

运行入口为 `run_ag607wt_case.py`；两套输出分别写入 `runs/full_1x/` 和 `runs/full_2x/`。`run_state.json` 指向最后的有效检查点，可用同一个 pair 脚本加 `--resume-existing` 恢复。控制器会把每网格约每 0.5 ps 及首次液相、连通熔池、稳定池和关热检查点推送到目标分支。2× 短测耗时约 32 秒/1 fs；完整 2× 结果需要长时间运行。

两套完整计算结束后，`compare_ag607wt_meshes.py` 会选择共同检查点时刻，对 2× 字段做 2×2×2 块平均后与 1× 比较，并输出 Ag 区温度误差、液态体积分数和最大连通液态团簇数据；`render_ag607wt_results.py` 会生成温度/液相截面和熔化—凝固历程图。报告会标记两网格处于不同加热/冷却阶段的时刻；单次网格配对不会被解释为收敛证明。

## 输入来源与旧算例状态

- 原 QSGS 父几何 SHA256：`7e0e1509b414597ce4e2a1c9fc387ed9656caf507f1186e2d28de79c3c1bfc43 `
- 本次新几何 SHA256：`edebee1d8d254a60e51c7de0b63d96c10263020989c9c963449b932323e95695`
- Ag70 2× 旧任务已停止；其 1.31 ps 检查点恢复验证通过并已推送。旧目录、几何、结果和检查点均保留。新配比的连续任务使用独立输出目录。

