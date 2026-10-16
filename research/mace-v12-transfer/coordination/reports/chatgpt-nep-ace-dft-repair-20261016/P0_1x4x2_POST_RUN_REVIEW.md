# P0：32原子首次实际DFT计算与历史标签比较审查

日期：2026-10-16。**本文件为ChatGPT本聊天界面读取现有GitHub原始数值并做离线对比后生成，未执行新DFT。**

## 数据来源与证据

- 新DFT结果分支：`cloud-h-p0-20261016`，提交 `ce501377a33a36b059621adefbf2c18ec9b6a2f6`。
- 结果目录：`research/mace-v12-transfer/coordination/reports/chatgpt-nep-ace-dft-repair-20261016/CLOUD_DFT_RESULTS_20261016/H_Ti_gap_32_k1x4x2/`。
- 新 `DFT_P0_RESULT.json` 原始SHA-256：`0a6fe5669af67200beea1af442d078408e277c850b5e1464b8e9313b062999df`（由云端归档提供）。
- 历史数据分支/提交：`499ddf0870235594444abad780582ea4fcc89d61`：
  `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v16_interface_energy/data/train.extxyz` 第0帧，以及
  `research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v17_rank_recovery_parallel/calculations/Ti_gap_2p3_rank_recovery_v17_01/`。
- 本次比较直接解析历史源extxyz及新32×3逐原子力；三个数据源32原子元素顺序和笛卡尔坐标逐值一致（历史extxyz所保存的小数精度内），没有重新映射/松弛。

## 数值对照（原子力均按行索引对应，单位eV/Å）

| 对照 | Native E（eV） | 相对新1×4×2总能差（eV） | 新减旧逐原子力向量RMS（eV/Å） | 最大原子力向量差（eV/Å） |
|---|---:|---:|---:|---:|
| 历史1×4×2（居中未验证） | −207.463968321659 | +0.600731579425 | 0.101692292030 | 0.249779232034 |
| 历史Γ单点（1×1×1） | −212.001510218357 | +5.138273476123 | 2.794235592585 | 5.761788170843 |
| **本轮Γ居中1×4×2（新DFT）** | **−206.863236742233** | **0** | **0** | **0** |

注意：上表“相对新1×4×2总能差”表示“新E − 历史E”，并非每原子值；和历史1×4×2之差为18.772861857 meV/atom。

### 本轮实际完成的DFT

- SCF 53迭代收敛；32×3力均有限；
- PBE/PW 500eV/Fermi0.10eV/Γ居中1×4×2，4 MPI进程；GPAW26.7.0、ASE3.29.0；
- Native E = **−206.8632367422333 eV**；
- Force-consistent free energy = **−207.10796845099478 eV**；
- Fmax = **3.3727980494371663 eV/Å**，冻结结构本身不是力平衡构型；这不影响固定几何SCF成功的事实；
- 用时914.774秒（约15分15秒）；
- 新旧文件内容的SHA可核验，但SCF收敛≠k网格收敛。

## 结论

```text
FIRST_CLOUD_SCF_CONVERGED=YES
NEW_DFT_PHYSICAL_LABEL_PRESERVED=YES
HISTORICAL_KPOINT_CENTERING_VERIFIED=NO
HISTORICAL_1x4x2_FORCE_AGREEMENT_WITH_NEW=NO
NEW_GAMMA_CENTERED_1x4x2_TO_1x6x4_COMPARISON_AVAILABLE=NO
REFERENCE_KPOINT_CONVERGENCE=NOT_VALIDATED
NEP_ACE_PRODUCTION_TRAINING_ALLOWED=NO
```

历史1×4×2与新1×4×2在名义网格相同但无法确认k点采样点完全相同；它们的能差不能直接称为k点收敛差。不要混用旧Gamma重算数据与当前1×4×2构建同口径能量训练集。

## 建议唯一下一算例：Γ居中1×6×4

1. 以新1×4×2的**完全相同原子位置、晶胞、PBC、PBE/PW500、Fermi0.10、PAW文件SHA、收敛条件、非自旋极化设置、4 MPI rank**为基线，**只改k点网格为Γ居中1×6×4**。
2. 检查实际BZ/IBZ点、收敛次数、墙钟时间、峰值内存、原子力(N×3)、native总能与force-consistent自由能；保存日志和SHA。
3. 仅比较这两次新DFT（同一k点居中定义）：`|ΔE_native|/32 ≤2 meV/atom`，`|ΔF_free|/32 ≤2 meV/atom`，`RMS_i ||F_i^(1x6x4)-F_i^(1x4x2)|| ≤0.01 eV/Å`；额外报告最大原子向量差。
4. 若三项皆通过，只能暂定本结构的k点差异在门槛内；仍需说明其他数值设置敏感性，不批准NEP/ACE生产。若有任一失败，先报告并停，**不要自动运行1×8×6**。
5. 不覆盖云端P0结果；将新标签归档到`cloud-h-p0-20261016`或新独立结果目录，完成后停止并返回commit/SHA供本窗口复核。

**本文件是下一步审核意见，不是已经发起的云端作业。**
