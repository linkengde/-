# 200 × 200 × 300 Å 模型：Ag 重新加热及熔化判定

## 结论

这次探索性计算支持：**Ag 在高温阶段呈液态/熔融态特征**。在 0.8 ps（Ag 温度约 1803 K）时，OVITO PTM 仅识别到 3 个 FCC Ag 原子、没有 HCP 原子；到 1.1 ps（约 1967 K）PTM 已不再识别到 FCC/HCP/BCC/ICO 晶体结构。与此同时，第一配位壳 Q6 和晶格位移 MSD 持续下降/上升。

这不能给出准确的熔化起点或熔点：初始 300 K 构型本身只有 7.76% Ag 被 PTM 识别为 FCC/HCP，其余多数为不规则局部环境。因此本结果确认的是**加热后的高温 Ag 已有强烈液态特征**，不是从完整晶态银测得的熔点。该模拟是对 Ag 直接均匀 NVT 加热，**没有使用 TTM**。

## 模型和计算设置

- 盒子：200 × 200 × 300 Å，458,248 个原子，其中 Ag 106,396 个。
- 数据起点：`data/ag_equilibrated_300K.data`，300 K 热化后的坐标和速度。
- 原子类型：Ti=1、Si=2、C=3、Ag=4。
- Ag–Ag：仓库中的 `models/ag_ti3sic2_qsgs_1to1_finer/potentials/Ag_u3.eam`；SHA-256 为 `cb3d605cd27ec6efb50c4dd7462c2ab3d967e63775dc4400de9134315a44d807`。
- 骨架：Ti、Si、C 原子坐标完全固定，作为几何支撑/障碍；它不参与动力学、吸热、散热或晶格传热。
- Ag–骨架：LJ/WCA 截断排斥项，ε=0.5 eV；Ti–Ag σ=2.08210927 Å、截断=2.33708864 Å；Si–Ag σ=2.14179034 Å、截断=2.40407837 Å；C–Ag σ=1.91771227 Å、截断=2.15255924 Å。LJ 势能在截断处平移为零。
- 边界：x/y 周期，z 非周期；底部 Ag 使用 `wall/reflect` 防止越过 z 下边界。
- 积分：LAMMPS `metal` 单位，时间步长 0.001 ps（1 fs），Ag 使用 Nose–Hoover NVT，阻尼时间 0.02 ps。
- 加热：0–0.5 ps 从 298.72 K 升至 1800 K；0.5–1.0 ps 保持 1800 K；1.0–1.2 ps 升至 2400 K；1.2–1.7 ps 保持 2400 K。
- 保存：Ag-only 轨迹每 100 步（0.1 ps）一帧，另存 1.7 ps 全体系快照。

以上设置是低成本探索模型。固定骨架及简化排斥交叉项不能代表真实 Ag/Ti₃SiC₂ 界面粘附、化学作用或骨架传热。

## 熔化判据和结果

Q6 用 LAMMPS `compute orientorder/atom degrees 1 6 cutoff 3.8` 在轨迹坐标上重算，截断半径 3.8 Å；表中 Q6>0.30 仅为描述性统计线，不单独作为熔化阈值。PTM 用 OVITO 3.17.2 对 Ag-only 轨迹逐帧分类。MSD 是 LAMMPS 对 Ag 组相对初始热化构型计算的累计 MSD。

| 时间 (ps) | Ag 温度 (K) | 平均 Q6 (3.8 Å) | Q6>0.30 比例 | PTM FCC+HCP 比例 | Ag MSD (Å²) |
|---:|---:|---:|---:|---:|---:|
| 0.0 | 298.7 | 0.2053 | 45.40% | 7.760% | 0.000 |
| 0.5 | 1743.2 | 0.1754 | 40.29% | 0.215% | 1.124 |
| 0.8 | 1802.9 | 0.1335 | 29.71% | 0.00282%（3 个 FCC 原子） | 2.914 |
| 1.1 | 1966.5 | 0.1213 | 26.58% | 0 | 4.575 |
| 1.7 | 2395.8 | 0.0920 | 18.79% | 0 | 8.952 |

末段 1.2–1.7 ps 的 MSD 增加约 3.77 Å²，且 PTM 全部落入 Other，支持 Ag 在该高温阶段处于强烈无序、可迁移的状态。Other 单独并不等于液体；此处是与 Q6 衰减、MSD 增长和高温共同判断。

计算结束仍有 458,248 个原子（Ag 106,396 个），没有丢原子。末态 Ag 最大瞬时力约 16.85 eV/Å；这项探索性结果未做力学平衡或能量守恒认证。

## 文件和复现

- `inputs/in.reheat.lmp`：重新加热的 LAMMPS 输入。
- `inputs/in.reanalyze_q6.lmp`：在相同坐标轨迹上用 3.8 Å 第一壳截断重算 Q6。
- `data/ag_equilibrated_300K.data`：复现加热需要的起始坐标和速度。
- `trajectories/ag_reheat_q6_ptm_msd.dump`：逐帧 Ag 轨迹；Git LFS 存储。
- `trajectories/ag_reheat_q6_ptm_msd_final.dump`：1.7 ps 全体系快照；Git LFS 存储。
- `logs/`：本次 LAMMPS 实际输出日志。
- `analysis/`：OVITO PTM 统计和结果合并脚本。
- `results/reheat_q6_ptm_msd_summary.csv`：逐帧合并的温度、Q6、PTM、MSD、最大力和原子数。

需要已安装 LAMMPS（含 EAM、LJ、orientorder/atom）和 OVITO Python 模块。示例：

```bash
LMP_BIN=/path/to/lmp ./run_reheat.sh
OVITO_PYTHON=/path/to/ovito-python ./run_analysis.sh /path/to/ag_reheat_q6_ptm_msd.dump /path/to/screen.reheat_q6_ptm_msd.out
```

脚本默认新建时间戳结果目录，不覆盖归档轨迹。软件版本及文件哈希见 `SOFTWARE_AND_HASHES.txt` 与 `SHA256SUMS`。

## 数据版本校验

`SHA256SUMS` 列出本归档内的输入结构、轨迹、日志、脚本和结果的 SHA-256；EAM 势文件的仓库相对路径及 SHA-256 见 `SOFTWARE_AND_HASHES.txt`。首次下载含 LFS 轨迹的分支后，执行 `git lfs pull` 可取回大文件内容。
