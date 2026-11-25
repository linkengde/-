# Ag/Ti₃SiC₂ TTM 喷溅探索：25.5 ps 检查点

本目录保存 25.5 ps 的安全检查点、该帧的原子轨迹、喷溅/表面形貌统计、OVITO 三维视图，以及续算脚本。数据作为**探索性结果**归档，未被并入 `main`。

## 结果摘要

- 体系共 308,510 个原子，其中 Ag 71,630 个；TTM 耦合 Ag 50,722 个。
- 时间步长 0.5 fs，检查点为 25.500 ps / step 51,000；z 向盒子扩至 400 Å。
- 25.5 ps 时，Ag 热耦合子集平均温度为 4,477.69 K；电子温度最大值为 9,080.69 K。
- 按 Ag–Ag 距离 4 Å 连通、并要求整个连通团簇高于初始表面 z=150 Å，识别到 425 个脱离团簇、共 4,130 个 Ag（总 Ag 的 5.77%）；另有 29 个孤立 Ag，其中 25 个当帧 vz>0。这是模型中的喷溅/外逸筛查，不能单独证明热蒸发。
- 中心 z90 比 20–40 Å 环带高 6.68 Å，显示中心隆起，未显示熔坑。
- 后续段在 25.9695 ps 触发保护停止，瞬时 Fmax=50.6254 eV/Å。未完成段没有新轨迹帧或可用重启；最后安全检查点仍为 25.5 ps。不要将部分 CSV 当作完成结果，也不要从其恢复。

## 模型与适用范围

此次探索使用纯 Ag Foiles EAM；Ti/Si/C 骨架原子完全固定；Ag–骨架相互作用使用 WCA 排斥；并采用 Python TTM 参考耦合。它没有经过验证的 Ag/Ti₃SiC₂ 四元素化学界面势。骨架形貌为当前 `thinner_scaffold_candidate_03` 版本，中心区域看起来偏团块；本归档保留原始几何，后续更细骨架应另建模型版本。

这些结果可用于复现本模型观察到的原子运动和外逸筛查，不应作为真实触头蒸发率、熔坑尺寸或定量材料预测。

## 归档文件

- `checkpoint_25p5ps/`：相互匹配的 LAMMPS restart、TTM 网格状态、25.5 ps dump、能量/温度指标及校验元数据。
- `failed_attempt_25p5_to_26p0ps/partial_continuation_metrics.csv`：保护停止前记录到 step 51920（25.96 ps）的部分指标；不是检查点。
- `analysis/`：进度报告、逐检查点统计 JSON、单帧分析程序、OVITO 渲染程序和 PNG 图。
- `analysis/q6_ptm_reference/`：OVITO 着色所需的 Q₆/PTM 分析辅助程序及几何参数。
- `reproduction/continue_heat_extend_z.py`：从成对 TTM/LAMMPS 检查点续算的脚本。
- 仓库模型目录中的 `potential_validation/geometry_constraint_B/thinner_scaffold_candidate_03_20261125/ag_ti3sic2_qsgs.data`：检查点元数据引用的原始结构输入，SHA-256 与原运行一致。

## 校验和

在本目录运行：

```bash
sha256sum -c SHA256SUMS.txt
```

## 复算 25.5 ps 单帧统计

需要 Python、NumPy、SciPy：

```bash
python analysis/analyze_checkpoint_frame.py
```

## 重新渲染 OVITO 图

需要与原运行兼容的 OVITO Python、PySide6、Vulkan 驱动（无 GPU 时可用 Mesa lavapipe）：

```bash
python analysis/render_splash_crater_25p5ps.py
```

输出写入 `analysis/`。渲染只读取归档的 25.5 ps dump。

## 从检查点继续的依赖关系

续算脚本使用仓库中的 `models/ag_ti3sic2_qsgs_1to1_finer/ag_ttm_reference.py`、Ag EAM、TTM 参数表、WCA 骨架边界输入和本次归档的初始结构文件。其 `MODEL_DIR` 常量默认指向 `/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer`；在其他 checkout 路径运行时应先按 README 所述改为本地仓库路径。

脚本默认执行 2,000 步，并会在 Fmax 超过 50 eV/Å 时保护停止。当前 25.5 ps 检查点的原加热设置在 25.9695 ps 已触发该门槛。归档不包含绕过保护阈值的设置，也不把失败段伪装成完整重启。继续研究前应先单独诊断异常高力局部环境；原始安全检查点保持不变。

`analysis/splash_crater_checkpoint_analysis.json` 保存 10.0–25.5 ps 检查点汇总。为控制仓库体积，本目录不包含全部历史分段的原子帧；只有 25.5 ps 原子帧随本归档上传。逐检查点汇总表不能仅凭本归档中的终点帧重新生成。
