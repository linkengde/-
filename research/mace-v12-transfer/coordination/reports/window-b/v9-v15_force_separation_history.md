# Ag–X 界面受力与中心分离力历史汇总（v9–v15）

**目的：** 给后续模型迭代整理可复用的 DFT 受力、中心接触对分离力和逐原子误差证据，明确哪些版次在哪些结构上未过筛。
**数据快照：** GitHub `main` 提交 `659c2f8dba6a7dbd9cce039c1db4f5c702d7250c`（2026-10-05 03:06 UTC）。
**范围限制：** 下文所有“误差”均来自仓库中已发表的 DFT 对照评估。V15 的当前仓库目录是残差定向 DFT acquisition；在上述快照中没有 v15 MACE checkpoint 或 v15 受力/分离力评估文件。若别的窗口已有未推送的 v15 模型结果，本汇总不包含它。

## 验收标准与指标含义

`coordination/ITERATION_PROTOCOL.md` 记录的暂定筛选门槛为：

- 能量绝对误差或每界面能量 MAE ≤ 10 meV/atom；
- 全结构原子受力**向量** RMSE ≤ 0.05 eV/Å；
- 标记中心 Ag–X 接触对的分离力绝对误差 ≤ 0.10 eV/Å。

中心分离力是沿标记原子对连线的相对受力分量；它与全结构的向量 RMSE 是不同指标。某个接触对的分离力通过，不能抵消周围 C/Ti/Si 原子较大的向量力误差。局部 RMSE 的定义是标记对任一原子 4.5 Å 内的原子，采用最小镜像距离。

## v13 初始独立界面筛查

v13 使用 v12 数据并增加 AgC 2.30 Å、AgSi 2.60 Å、AgTi 2.60 Å 三个固定距离 PW-PBE 点；v12 的 valid/test 未改动。对 v13 的独立 acquisition 结构按界面汇总如下。AgC、AgSi 各有两个结构，AgTi 有一个结构；表中的分离力误差是该界面结构中的最大绝对误差。

| 界面 | 最大能量误差 (meV/atom) | 力向量 RMSE (eV/Å) | 最大分离力误差 (eV/Å) | 筛选结论 |
|---|---:|---:|---:|---|
| Ag–C | 4.877 | 0.11368 | 0.35750 | FAIL：力、分离力超门槛 |
| Ag–Si | 2.537 | 0.09037 | 0.16178 | FAIL：力、分离力超门槛 |
| Ag–Ti | 4.561 | 0.06259 | 0.04392 | FAIL：力超门槛 |

这解释了为什么 v13 即使在相关的中间距离插值点上看起来改善，仍被标为 `FAIL_SCREENING`：更换 registry/局部构型后，Ag–C 和 Ag–Si 的受力/分离响应没有达到门槛，Ag–Ti 的整结构力 RMSE 也超过 0.05 eV/Å。

## v14 保留盲留出上的逐界面比较

下表使用同一组三条标记为 v14 holdout 的结构，回算比较 v9–v14。它是跨版本的同构型比较，不应倒称为 v9–v11 当年的盲测。另有数据完整性审计发现，Ag–Ti 这条所谓 v14 holdout 与 v14 训练帧 `AgTi_2p7A_LCAO` 在原子 ID、元素和坐标上完全相同；因此 Ag–Ti 结果只能作为该特定几何的参考标签一致性/模型评分，不能作为独立几何泛化证据。Ag–C、Ag–Si 仍是 v14 指定留出；其是否曾进入 v9–v11 的训练集，本报告没有逐版重建，故也不声称对旧版本盲测。每格格式为“力向量 RMSE / 中心分离力绝对误差”，单位均为 eV/Å。

| 模型 | Ag–C 留出 | Ag–Si 留出 | Ag–Ti 留出 |
|---|---:|---:|---:|
| v9 | 0.1843 / 0.2717 | 0.1445 / 0.1048 | 0.1966 / 0.1749 |
| v10 | 0.1880 / 0.1701 | 0.1468 / 0.0961 | 0.1933 / 0.2715 |
| v11 | 0.1870 / 0.1897 | 0.1393 / 0.0911 | 0.1814 / 0.3073 |
| v12 | 0.2126 / 0.0294 | 0.1366 / 0.3994 | 0.1414 / 0.1972 |
| v13 | 0.1633 / 0.0300 | 0.1189 / 0.2565 | 0.1681 / 0.2584 |
| v14 | 0.1537 / 0.0614 | 0.1178 / 0.2349 | 0.1466 / 0.1926 |

另附机器可读的 [V9–V11 误差明细 CSV](v9_v10_v11_force_separation_history.csv)：共 11 个结构、每个结构各含 V9/V10/V11 一行。它从已发布的 V12 候选比较表提取 8 个既有 gap、距离扫描和 midpoint 结构，并加入 3 个 V14 冻结留出上对 V9–V11 的回算结果。能量统一记录绝对误差；无标记 Ag–X 对的 gap 结构，分离力误差留空。CSV 每行保留源结果文件路径及能量、全结构力、最大力误差、分离力误差和门槛判定。

请注意：V14 冻结留出上的 V9–V11 数值是对相同新结构的回溯评分，便于横向比较；不能称作 V9–V11 当年独立盲测。部分旧 midpoint/距离扫描也来自相近母构型，不能替代全新 registry 验证。

解读：

- 这三条 v14 指定评分记录上，六个版本的全结构力 RMSE 全都高于 0.05 eV/Å；因此不能仅凭能量误差改善认定受力通过。Ag–Ti 记录存在上段所述的训练/测试几何重复。
- v14 的 Ag–C 中心分离力误差 0.0614 eV/Å 通过 0.10 门槛，但其全结构力 RMSE 仍为 0.1537，4.5 Å 局部 RMSE 为 0.1762，最大单原子误差 0.3241 eV/Å。问题主要在局部壳层原子，不只是中心 Ag–C 这一对。
- v14 的 Ag–Si 分离力误差 0.2349、Ag–Ti 为 0.1926，均未通过；三个界面的全结构力 RMSE 也都未通过。
- v14 能量在三个界面留出上均低于 10 meV/atom，但冻结测试集仍有 Si-gap 能量误差 15.91 meV/atom，以及 C-gap 力 RMSE 0.08263 eV/Å，故总判定 `SCREEN_FAIL`。

v14 留出上的局部误差定位：

| 界面 | 全结构 RMSE | 4.5 Å 局部 RMSE | 区域外 RMSE | 最大单原子力误差 | 中心分离力误差 |
|---|---:|---:|---:|---:|---:|
| Ag–C | 0.1537 | 0.1762 | 0.0937 | 0.3241 | 0.0614 |
| Ag–Si | 0.1178 | 0.1266 | 0.1079 | 0.2793 | 0.2349 |
| Ag–Ti | 0.1466 | 0.1775 | 0.0834 | 0.4652 | 0.1926 |

最大残差原子也不总是标记的 Ag–X 原子。例如 v14 Ag–C 留出中，最大力误差出现在 Ti#604590（0.3241 eV/Å），其后是 Ti#604806（0.2978）和 C#604602（0.2875）；中心 Ag–C 分离力本身却达到阈值。Ag–Ti 留出最大误差为邻近 C#480805 的 0.4652 eV/Å。这就是残差采样需要覆盖接触第一壳层和骨架邻居的原因。

仓库里的 v15 修正计划记录，v14 训练到 optimizer epoch 79，选择 epoch 71；训练配置已设 `forces_weight=1000`、`energy_weight=100`。因此，单独继续提高力权重不是有证据支持的首要修复。计划优先补不同 Ag registry 的独立留出和接触壳层/邻近骨架 DFT 标签，再用固定训练配置比较数据扩充效果。

## v12–v13 的相关中间距检查不能代替新构型验证

v13 在 Ag–C/Ag–Si midpoint 上相对 v12 的力 RMSE 明显下降，但 Ag–C 分离力误差仍约 0.110 eV/Å、Ag–Si 约 0.143 eV/Å。两个 midpoint 与后来加入训练的 AgC 2.30 Å、AgSi 2.60 Å 点来自同一母构型，只改变中心距离，属于相关插值检查，不是独立 registry 验证。

后续 v13 lateral-registry 检查显示：Ag–C 力 RMSE 0.08913、局部 RMSE 0.10147、分离力误差 0.17487；Ag–Si 对应 0.07689、0.07947、0.12866；Ag–Ti 的力 RMSE 0.06259、分离力误差 0.04392。它们说明增加距离点改善了部分插值指标，但尚未让所有接口同时满足向量力门槛。

## V15 当前状态与数据独立性

在本报告所用仓库快照中，V15 还不是已训练和已评估的势函数版本，而是三条残差定向 PW-PBE 采集任务：

- `AgC_residual_shell_v15_01`：60 轮 SCF 收敛并通过归档检查；输入几何中心 Ag–C 距离 2.2571 Å。摘要里的 DFT 分离力 −2.6244 eV/Å 是**参考值，不是 MACE 误差**。
- `AgSi_residual_shell_v15_01`：快照时正在计算，记录到 SCF 第 10 轮；
- `AgTi_residual_shell_v15_01`：快照时排队。
- 该快照没有 `mace_periodic_v15...` checkpoint 或 v15 force/separating-error assessment。因此“V15 模型未达标”的具体数值需从产出这些指标的后续提交补入，不能从 DFT acquisition 的原始力值推算。

已发布的 holdout-design audit 还发现，三个 V15 输入都是 v14 训练帧的局部小扰动：Ag–C 位移 RMS/max 为 0.0288/0.0602 Å；Ag–Si 为 0.0520/0.1270 Å；Ag–Ti 为 0.0440/0.0886 Å。它们可用于补训练残差，但不能作为 v15 盲测。审计在当前记录中没有找到三种界面可用的全新 registry/独立来源构型。

## 对后续训练的建议

1. 保留 v12、v13、v14 数据、模型和盲测原样，不覆盖旧结果。
2. 等三条 V15 DFT 标签全部通过哈希、结构、收敛及有限值校验后，才构建 v15 训练数据；明确哪些 V15 acquisition 会进入训练。
3. 在训练 v15 前冻结并哈希新的 blind holdouts。必须使用不同 registry 或独立来源的局部界面几何；同一母结构上的很小扰动不可算作独立验证。
4. 对新模型逐结构同时报告能量误差、全结构力向量 RMSE、接触对分离力误差、最大单原子误差及局部/区域外 RMSE。不要只报告全体系平均或一个中心对力。
5. 若新构型的 Ag–C/Ag–Si 仍失败，按误差最大原子及邻居环境补 DFT；不通过前不做候选大模型长 MD/TTM。

## 可复核证据

- 门槛和数据流程：`coordination/ITERATION_PROTOCOL.md`
- v13 独立 acquisition 评分：`periodic_interface_v4/mace_periodic_v13_interface_energy/results/v13_independent_holdout_comparison.json`
- v13/v14 原子级局部误差：`coordination/reports/window-b/v13_force_localization.{md,json}`、`v14_force_localization.{md,json}`
- v14 筛查：`periodic_interface_v4/mace_periodic_v14_interface_energy/results/v14_validation_assessment.json` 和 `v14_independent_holdout_comparison.json`
- V14 Ag–Ti 完全重复及 LCAO/PW-PBE 受力不一致的限制：`periodic_interface_v4/V15_FORCE_SEPARATION_PLAN.md`（“Dataset-integrity finding”）和 `WORK_LOG_2026-10.md`
- 本窗口整理的 V9–V11 机器可读明细：`coordination/reports/window-b/v9_v10_v11_force_separation_history.csv`
- V15 力/分离力修正计划：`periodic_interface_v4/V15_FORCE_SEPARATION_PLAN.md`
- V15 采集进度：`coordination/progress/window-a.json`、`periodic_interface_v4/pbe_interface_v15_targeted_acquisition/partial_archive_manifest.json`
- V15 几何与 holdout 独立性：`coordination/reports/window-b/v15_holdout_design_audit.md`
