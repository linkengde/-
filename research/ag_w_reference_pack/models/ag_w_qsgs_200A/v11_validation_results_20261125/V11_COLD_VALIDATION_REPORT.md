# Ag–W V11 验证与 TTM 软件联调记录

日期：2026-11-25（本地验证）
分支：`research/ag-w-ideal-gyroid-v11-20261124`
生成几何提交：`e97c652d95cc9ec52c469a588cb748197396325f`
静态/冷态输入来源：`84467c5e`
范围：仅新增本验证目录和 V11 TTM 求解器/工作流更新；没有修改 `v11_generated_data/` 或 V5–V10。

## 结论

| 项目 | 状态 | 说明 |
|---|---|---|
| Gyroid 几何生成 | PASS | 462,180 原子，Ag/W 质量比 0.999995，双方连续 |
| Dai EFS 静态筛查 | 探索性 PASS | FIRE 22 步后以能量容差停止；残余力仍高 |
| 300 K 冷态基线 | **FAIL** | Ag/W 温差 82.9 K；再续跑 1 ps 后扩大到 92.5 K，势能仍下降 |
| 10 Å 真实几何映射 | PASS（仅映射预检） | 联合填充分数不超过 1；没有执行 TTM 步 |
| 无热源 TTM 联调 | PASS（软件联调） | 0.08 ps，联合能量账和网格能量账通过 |
| Ag 高斯热源 TTM 联调 | PASS（软件联调） | 累计热源 66.105 eV，能量账通过 |
| Ag 熔化/凝固物理结果 | **未获得** | 本轮没有验证晶体熔化、液相形成或凝固；不应称为熔焊模拟结果 |

因此，V11 的映射和 TTM 软件链路已能在这份真实几何上运行，但冷态基线仍不合格。本轮的 TTM 结果是非平衡起始态上的短时软件联调，不是定量材料预测。

## 几何、静态和冷态结果

结构有 Ag 291,274 个、W 170,906 个，总计 462,180 个；Ag/W 质量比 0.999995。Dai EFS FIRE 最小化在 22 步后因能量容差停止，最大移动原子力为 2.636 eV/Å，RMS 为 0.158 eV/Å，所以只记为探索性静态筛查通过，不代表低残余力平衡。

首段冷态过程：固定 z=0–20 Å，其他原子可动；300 K NVT 2 ps，再做 0.05 ps NVE。NVT 末端整体为 301.976 K，但 Ag/W 分别为 271.602/353.583 K，温差 81.981 K。NVT 最后 0.2 ps 势能仍下降约 250 eV。NVE 末端温差 82.860 K，日志精度内能量漂移为 0，原子数不变。

从 step 8200 诊断 restart 续跑 1 ps 共用 NVT + 0.05 ps NVE 后，NVT 整体温度 298.813 K，Ag/W 为 265.269/355.994 K，温差 90.725 K；续跑末端势能在最后 0.2 ps 仍下降约 165 eV。NVE 末端温差为 92.508 K。延长 NVT 没有消除相间温差。

LAMMPS 曾报告 `Temperature for fix modify is not for group all`。积分虽完成，后续仍需复核温度 compute 与 NVT 组设置。温差原因尚未唯一确定，不能直接归因于界面热导或势函数。

## 10 Å 映射预检

映射使用提交的 Gyroid 分数场和 step 8200 冷态坐标，命令使用 `--mapping-only`，没有原子积分、TTM 步或热源。

- 状态：`AGW_MAPPING_PREFLIGHT_PASS`
- Ag+W 最大几何填充分数：1.0；相分数联合容量未超限。
- Ag/W 原子到网格权重行和误差：2.22×10⁻¹⁶。
- Newton 层候选原子 367,582 个，其中 Ag 231,274、W 136,103；z<40 Å 的 205 个 Ag 被排除，z>200 Å 为 0。
- Ag 源空间积分功率：6,210.164 eV/ps；映射前后功率误差为 0（浮点精度内）。
- `run 0` 力筛查最大力 5.836 eV/Å；报告明确标注这不是平衡验收。

## 无热源与高斯源 TTM 软件联调

两组都从同一个 step 8200 restart 启动，使用同一个 10 Å 映射、Dai EFS、相同时间参数和随机种子；分别将热源倍率设为 0 和 1。电子初温按 Newton 层 Ag/W 各自去除质心运动后的平均晶格温度初始化，而不是把非平衡 restart 强设为均匀 300 K。初始 Ag/W 电子温度分别为 271.836/357.539 K。

设置为 dt=0.0001 ps（0.1 fs）、总耦合时间 0.08 ps；高斯 FWHM=0.01 ps、中心=0.02 ps，脉冲后冷却 0.02 ps；空间热源沿用 I₀=5 eV/(ps·Å²)、R=30 Å、穿透深度 20 Å、中心 (100,100) Å、表面 z=200 Å。Ag/W 电子界面热导采用临时探索值 1 GW/(m²·K)。

| 指标 | 无热源 | Ag 高斯热源 |
|---|---:|---:|
| 联调状态 | PASS | PASS |
| 累计源能量 | 0 eV | 66.105 eV |
| 最大联合能量残差 | 0.00727 eV | 0.00727 eV |
| 联合账容差 | 0.03197 eV | 0.09807 eV |
| 恒温浴能量 tally（正值表示能量进入恒温浴） | −31.970 eV | −31.969 eV |
| 最大网格步能量残差 | 1.43×10⁻¹⁴ eV | 7.33×10⁻¹⁴ eV |
| Ag 最大电子温度 | 311.5 K | 630.8 K（20.8 fs）|
| W 最大电子温度 | 357.6 K | 370.9 K |

热源组的累计输入为 66.105 eV。两组恒温浴能量 tally 均约 −31.97 eV；按 LAMMPS tally 符号约定，这表示恒温浴在该短跑中净向原子系统输入约 32 eV，与激光输入同量级。因此 source/no-source 对照可用于软件联调，但这段运行并非只有激光源在改变系统能量。报告出的峰值是**电子温度**，本轮没有记录或验证液相原子结构，也没有分析局部有序度；不能由此声称 Ag 已熔化。热源相分数和原子到网格映射在整个短跑中保持固定，无法描述大幅熔池扩展或相界移动。界面电子热导仍是未经 Ag/W 实测验证的假设值。

### 能量账实现诊断

本机 Python-LAMMPS 绑定最初找不到 `libmpi.so.12`；将匹配的 MPICH 库目录加入 `LD_LIBRARY_PATH` 后，Python-LAMMPS 接口启动正常。TTM 脚本在 `thermo 0` 时用 `get_thermo("pe")` 会因能量未计入而失败；V11 求解器现通过 `pair_modify compute yes` 和直接 PE compute 读取能量，避免打开会导致本机 Python/MPICH build 崩溃的 thermo 输出路径。

另外用 100 个无相互作用粒子独立核对了 Langevin `tally yes` 符号：粒子动能增加 11.55165 eV 时 fix tally 为 −11.55109 eV，说明 fix tally 与原子系统能量变化方向相反。TTM 联合能量账按该符号记账，并在无源/有源两组中闭合至上述残差。核对输入与日志见 `v11_langevin_tally_sign.in/.log`。

## 数据来源与哈希

V11 压缩结构 SHA-256：`ae09e6fa39947d3988d5af943af14462e0e10d330aa4941ab2905b287fb41b44`。解压后的原子数据 SHA-256：`475691478b020c8dff3d8db28bb57ced8e8ae69c3c6ff9351428deb95a64354a`。几何分数 NPZ SHA-256：`d2806712257feb3c4e0ab8136de6864979d6f09410397fe25f86c5576be95018`。本地生成结构解压后的字节与已提交文件一致；gzip 元数据可造成压缩文件自身 SHA 不同。

## 归档文件

- `v11_static_metrics.json`、`v11_cold_metrics.json`：静态和冷态指标。
- `v11_cold_followup_metrics.json`、`v11_cold_followup_timeseries.dat`、`in.v11_cold_followup`：冷态续跑记录。
- `v11_mapping_preflight.json/.lammps.log`：真实结构 mapping-only 结果。
- `v11_ttm_no_source.json/.lammps.log`、`v11_ttm_gaussian_source.json/.lammps.log`：匹配的 TTM 软件联调结果及完整温度/能量历史。
- `v11_ttm_exploratory_manifest.json`：校验输入结构、分数场和 restart 的联调 manifest。
- `v11_cold.restart.gz`、`v11_cold.dump.gz`：冷态失败诊断 restart/dump，仅可用于联调，不是已通过的激光起始态。
- `v11_validation_manifest.json`、`SHA256SUMS.txt`：综合状态和文件完整性校验。
- `v11_ttm_pair_commands.md`：复现实验命令。

## 下一步

先解决 Ag/W 冷态相间温差、持续势能下降和高残余力问题；再用通过的冷态状态重新执行配对 TTM。真正研究熔化/冷却还需要明确足够的输入能量，并增加原子结构/液相识别和动态相界映射。不要覆盖 `v11_generated_data/`，也不要改写 V5–V10 历史数据。
