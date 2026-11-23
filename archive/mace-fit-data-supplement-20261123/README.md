# Ag/Ti₃SiC₂ 势函数数据补充（2026-11-23）

本归档接在已上传的探索性拟合归档之后：

- 基础拟合分支：`codex/mace-fit-ace-exploratory-20261119`
- 基础提交：`ad45eff650886243898bc8b0d3bee3402237603a`
- 基础目录：`archive/mace-fit-ace-exploratory-20261119/`
- 基础归档清单共 135 项，远端副本 SHA-256 全部通过。
- 基础 `SHA256SUMS.txt` 的 SHA-256：`fcb1d26411dadf12b6a130f788af6410b08f324cde190b8448b59363730724d4`。

## 本次数据补充

`ag_solid_controls/Ag_solid_controls_train_candidate_unassigned.extxyz` 汇集仓库既有的 7 个已收敛纯 Ag 固态 DFT 控制，覆盖基准、±0.01 Å 位移、±2% 体积和 ±0.5% 剪切。每个来源输出目录的 SHA-256、GPAW 日志及验证报告均已复核。它们来自同一个 Ag32/COD 9009647 结构族，只可作为一组训练候选，不能拆到验证/测试，也不代表有限温振动。

`ag_liquid_first_label/source/` 保存一个 32 原子纯 Ag 液态候选的生成输入、采样轨迹、结构筛查和 GPAW 计算脚本。Williams–Mishin–Hamilton 2006 Ag EAM 仅用于产生候选构型；它不作为能量或力参考标签。

DFT 使用 GPAW/PBE/PW500 eV、3×3×3 k 点、0.10 eV Fermi–Dirac 展宽和四个 MPI rank。候选结构来自 1500 K EAM NVT 轨迹的第 44000 步。训练数据整合保持关闭；只有 SCF、能量/力映射和独立标签审查都通过后，才可将该帧加入训练候选集。

## 当前结论

本归档不宣称已得到生产级势函数。新增标签属于纯 Ag 高温液态的探索性候选，不能替代固态 Ag、固液共存和高温 Ag/Ti₃SiC₂ 界面数据，也不能单独作为独立验证集。最终 DFT 状态、误差统计和文件哈希见 `ag_liquid_first_label/`。


## 2026-11-23 探索性 ACE 重拟合和冷态筛查

新增目录 `exploratory_ace_refit_20261123/` 保存基于 71/7 基础划分和 12 个新标签候选的 83/7 探索性 ACE 重拟合、五种正则化候选及冷态界面筛查。开发集来源族相关，不是独立测试集。虽然多组候选的开发误差满足宽松探索指标，五个新模型都没有通过 52 原子 Ag/Ti₃SiC₂ 冷态界面的松弛门槛，因此没有进入 NVE 或加热。旧 ACE_energy100 的冷态测试仅作为有限参考，不能证明它描述了液态 Ag。

`ag_liquid_density_reference_attempts/.../dft_attempt_002/` 保存密度修正 64 原子液态 Ag 构型的第二次 GPAW 尝试。它在一次未收敛 SCF 后主动中断，无检查点、无正式力标签，不加入拟合。完整范围和判定见探索性重拟合报告。
