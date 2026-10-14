# Stage E3 最终决策

## 状态

| 项目 | 状态 | 依据/限制 |
|---|---|---|
| E2_PACKAGE_SHA | PASS | 13 个文件存在；12 个载荷 SHA 全通过 |
| BASELINE_REPRODUCED | PASS（归档复核） | E2 已记录坐标匹配和 48/114 原子基线 RMSE；本轮未重算 |
| PARAMETER_MAPPING | PARTIAL | 已核验原文件 27×14 字段及数值共享模式；原论文作者约束未在本地验证 |
| PAPER_SOURCE_E_LOCAL | NOT_VERIFIED | 原文/附录未在仓库或附件中找到；仅接收外部 27 条/378 数值的核验结论 |
| PER_ATOM_PARAMETER_RESPONSE | NO | E2 只有汇总 RMSE，无参数变体逐原子力向量 |
| PARAMETER_IDENTIFIABILITY | NOT ESTABLISHED | 无逐原子响应 Jacobian；不报告秩或奇异值 |
| DATA_SUFFICIENCY | INSUFFICIENT | 缺体积/剪切 DFT stress、Si 位移、独立 clean surface、独立 DEV/TEST；平衡结构也未确认 |
| INDEPENDENT_VALIDATION_AVAILABLE | NO | 角色文件冻结 TRAIN/DEV/TEST 为空 |
| START_TRIAL_FIT | NO | 先补受约束逐原子 run 0 响应；补关键 DFT；由 C 冻结 split 后重新审核 |
| PRODUCTION_READY | NO | 无候选、无独立验证，也未验证应用温度范围 |

## 对问题的回答

1. **最少需要拟合哪些独立参数？** 现有证据不能选定某个具体行参数。第一轮最多开放一个已核验 tied block 中的一个标量；候选字段家族是 `A/B/lambda1/lambda2`，pair 候选是 Ti–C、Ti–Si、Si–C。当前不能把四个字段或三个 pair 一次放开。需要先用参数块逐原子力响应筛选确定成员。
2. **现有数据能否约束这些参数？** 不能达到行级可辨识性。现有数据能说明全列方向改变会影响聚合力 RMSE，且残差很大；它不能区分参数块、补偿关系或唯一根因。
3. **能否做探索性初步试拟合？** 当前不建议，也不启动。先执行经批准的逐原子 run 0 判别测试、补充必要 DFT 并冻结 family-level TRAIN/DEV/TEST。否则连探索方向都无法可靠选择。
4. **哪些额外 DFT 不可避免？** 为覆盖当前指定目标，至少需要体积±2%与剪切±0.5%共4点、Si 位移±0.02 Å共2点、独立 clean-surface 两种终止的参考/局部响应共4点，以及独立 DEV/TEST 结构家族。若缺可信平衡晶格/内部坐标，还需一次 PBE 体胞弛豫工作流或审核通过的等价来源。
5. **能否减少新增数量？** 可以把 DEV/TEST 各缩成一个独立结构点，将推荐14个静态点降到12个，但失去留出家族内响应检验。不要删除体积/剪切、Si、clean surface 项目，同时仍声称能覆盖对应目标。若放弃某项目目标则可以更少，但适用范围必须收窄。
6. **预计多少静态计算？** 推荐14个新的固定几何 DFT 单点；低配验证方案12个；另有条件性的一次体胞弛豫。6个 bulk 体积/剪切/Si 点按现有48原子 k10×10×4 单点的 11722.3 s、4 MPI 历史值粗估约19.5小时连续墙钟（约78 MPI-core-hours）；这是外推，不是保证。表面和独立验证模型时长缺少同规模参照，不能可靠估算总墙钟。
7. **何时可生成候选 Tersoff 文件？** 在原文件只读保留的前提下，来源/单位/共享约束核实、逐原子响应可辨识、DFT 数据和角色按家族冻结、有限参数候选在 TRAIN 上满足力/能量/弹性/晶格/表面目标、DEV 与 TEST 独立通过且无明显回归后，才可生成带新 SHA 的候选副本。候选仍须人工审查；未通过独立测试不得标为生产可用。

## 最小阻碍与后续顺序

最小阻碍是缺少逐原子参数块响应和独立验证角色。建议先单独授权 `PARAMETER_RESPONSE_LIMITS_AND_RUN0_PLAN.md` 的静态 `run 0` 屏查（不改变角色或势文件），并获取论文原文/附录以确认作者参数约束；然后由 C 审核结构角色和最小 DFT 清单。当前不优化参数，不启动 DFT/MD/TTM，不修改原始势函数。
