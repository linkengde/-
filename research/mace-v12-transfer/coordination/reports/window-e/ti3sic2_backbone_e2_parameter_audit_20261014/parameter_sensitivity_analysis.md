# TSCBOP 参数族静态敏感性

## 方法

以归档原参数文件 `TSCBOP.tersoff` SHA-256 `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7` 为只读源。使用 LAMMPS 22 Jul 2025 Update 4 对 48 原子基线 + 同一 parent 下四个 Ti/C 位移点、以及一个114原子界面派生 TSC-only slab 做 `run 0`。每个输入 0 个时间步，无积分、无坐标变化。MPI=1、OMP=1。基线力ID/元素和周期等价坐标匹配误差≤7.4×10⁻¹⁵ Å；基线RMSE复现既有48原子与114原子结果。

对Tersoff文件的13个连续数值列（`gamma, lambda3, c, d, h, n, beta, lambda2, B, R, D, lambda1, A`）分别生成临时参数副本，把该列所有非零三体条目统一乘以0.99或1.01；整数离散参数 `m` 固定。共26个副本、156个结构/副本 `run 0`，另有6个原参数基线。原文件SHA最终仍一致。该方向保留了同列原有相等值组、符号和通用正值范围；但原论文特有的跨列耦合约束没有从论文/来源文档中核验，所以这些是“全列方向”诊断，不能当作合格的行级独立参数化或候选参数集。

表中 `dRMSE/dlnp=(RMSE(+1%)-RMSE(-1%))/0.02`，单位eV/Å。负号表示这个全列方向在局部增大可降低此数据组RMSE；它不是因果证明或优化结果。

## 基线复核

| 结构族 | Ti RMSE | Si RMSE | C RMSE | 全部 RMSE (eV/Å) |
|---|---:|---:|---:|---:|
| 48原子基线及4个同父位移点，10×10×4 | 2.2206 | 0.0062 | 2.1290 | 1.9941 |
| 114原子单个界面派生TSC-only slab | 2.6242 | 2.6072 | 2.4433 | 2.5624 |

48原子家族的Si基线/小扰动RMSE很低，反映此结构和这些z位移下Si受力很小；它不能代表Si受扰动或表面环境。

## 最强局部方向

| 结构族 | 元素 | 参数列 | 基线 | -1% | +1% | dRMSE/dlnp (eV/Å) |
|---|---|---|---:|---:|---:|---:|
| bulk_parent_family | Ti | `lambda1` | 2.2206 | 2.8275 | 1.6152 | -60.615 |
| bulk_parent_family | Ti | `lambda2` | 2.2206 | 1.6796 | 2.7518 | 53.610 |
| bulk_parent_family | Ti | `A` | 2.2206 | 1.9659 | 2.4752 | 25.464 |
| bulk_parent_family | Ti | `B` | 2.2206 | 2.4533 | 1.9878 | -23.279 |
| bulk_parent_family | C | `lambda2` | 2.1290 | 2.1799 | 2.0811 | -4.944 |
| bulk_parent_family | C | `lambda1` | 2.1290 | 2.1533 | 2.1035 | -2.493 |
| bulk_parent_family | C | `A` | 2.1290 | 2.1043 | 2.1537 | 2.470 |
| bulk_parent_family | C | `B` | 2.1290 | 2.1333 | 2.1247 | -0.429 |
| interface_derived_slab | Ti | `lambda1` | 2.6242 | 3.3650 | 2.0305 | -66.726 |
| interface_derived_slab | Ti | `lambda2` | 2.6242 | 2.1425 | 3.1611 | 50.931 |
| interface_derived_slab | Ti | `A` | 2.6242 | 2.3841 | 2.8745 | 24.516 |
| interface_derived_slab | Ti | `B` | 2.6242 | 2.8518 | 2.4067 | -22.256 |
| interface_derived_slab | Si | `lambda1` | 2.6072 | 2.2861 | 2.9868 | 35.033 |
| interface_derived_slab | Si | `lambda2` | 2.6072 | 2.8808 | 2.3583 | -26.125 |
| interface_derived_slab | Si | `B` | 2.6072 | 2.5448 | 2.6716 | 6.337 |
| interface_derived_slab | Si | `A` | 2.6072 | 2.6518 | 2.5647 | -4.356 |
| interface_derived_slab | C | `lambda2` | 2.4433 | 2.5137 | 2.3778 | -6.796 |
| interface_derived_slab | C | `lambda1` | 2.4433 | 2.4766 | 2.4083 | -3.419 |
| interface_derived_slab | C | `A` | 2.4433 | 2.4077 | 2.4789 | 3.562 |
| interface_derived_slab | C | `B` | 2.4433 | 2.4550 | 2.4316 | -1.172 |

体相 Ti/C 的最大方向影响来自 `lambda1`, `lambda2`, `A`, `B`；`d,n,c` 次之。该结果只说明这些参数列全局缩放对当前force residual较敏感。两个元素并不总偏好同一方向：例如体相 `lambda2` 减小改善Ti、增大改善C；表面派生slab的Si对 `lambda1` 和 `lambda2` 的偏好又与Ti/C部分相反。单一全局因子不能同时修复所有组。`R,D` 在48原子基线/扰动族的±1%探针中不改变RMSE，说明这些几何未有效采样切换区；不代表截断参数正确。

## 根因结论

在固定基线结构上，原TSCBOP相对PBE Ti/C力残差为约2 eV/Å，远大于已归档DFT数值网格差异；在界面派生slab的Ti、Si、C残差均约2.4–2.6 eV/Å。此前分项结果显示，这主要是骨架内势误差，不是仅调Ag交叉项可以修复。参数族敏感性补充了“哪些方向影响残差”，但没有辨识单个行参数、唯一根因或参数可辨识性。

**因此：** 优先关注Tersoff径向强度/衰减族 `A,B,lambda1,lambda2`，不应同时自由优化；对角向/键级族 `c,d,h,beta,n,gamma,lambda3` 先保持冻结，等径向约束验证后再单独检查。此顺序只是小型诊断排序，正式参数块必须由原论文/实现约束确定。没有独立DEV/TEST，不能基于这些结果宣称泛化。
