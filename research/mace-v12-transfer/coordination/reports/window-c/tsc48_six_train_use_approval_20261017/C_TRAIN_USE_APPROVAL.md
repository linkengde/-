# C：六点有限探索性 TRAIN 用途批准

日期：2026-10-17（Asia/Shanghai）。申请：E commit `5a2480dad7b37ac707bf435a9201b28a1554abc9`。

REVIEW_STATUS = APPROVED

TRAIN_USE = APPROVED_LIMITED_EXPLORATORY_FORCE_ONLY

TRAIN_GEOMETRIES = 6；SOURCE_FAMILIES = 1；DEV = []；TEST = []。

FIT_AUTHORIZED = NO；EXECUTION_AUTHORIZED = NO。

## 批准范围

批准在新增、独立的用途映射中将以下六项登记为 TRAIN 资格，仅限后续受约束、力重点的Ti₃SiC₂骨架TSCBOP探索性修复准备：

- TSC48_VOL_m02
- TSC48_VOL_p02
- TSC48_BASAL_XY_SHEAR_m005
- TSC48_BASAL_XY_SHEAR_p005
- TSC48_SI_ID7_Z_m020
- TSC48_SI_ID7_Z_p020

逐项输入、RESULT及labelled.extxyz SHA、来源路径和批准范围见TRAIN_USE_APPROVAL.json。本映射记录当前批准用途；原始DIAGNOSTIC_ONLY档案、旧角色冻结文件、历史C审核和标签字节均不变。本轮没有复制标签或建立物理训练数据文件。

## 审核依据

已从真实GitHub fetch申请固定提交，申请两项清单SHA通过。申请中的六项精确身份与已审核C v2 manifest一致；本次再次读取申请提交中的六个输入/结果/标签文件，三类SHA逐项一致。申请所引用C报告、裁定、草案的文件哈希在其固定提交上复核通过，scope未扩大。

六点SCF和完整力标签已有C审查；执行锁、runner及批次PAW证据已有来源补证审核。现有证据足以接受这次明确限定的探索用途；批准不将未验证数值误差宣布为收敛。

## 强制用途边界

只使用原子力目标，单位eV/Å。energy和free_energy、stress保持诊断，禁止进入本次批准的初始loss，禁止据此标定或认证弹性常数。MP6×6×2扰动结构网格/应力未独立收敛、PAW批次级而非逐作业绑定、五项launcher退出码未知等限制随用途映射保留。

六项属于同一bulk_parent_48父族；六帧不等于六个独立家族。所有同族后代不得作为独立DEV/TEST。旧五点、H候选及其他后代不自动取得TRAIN用途。独立DEV尚缺；TEST由独立evaluator管理，本轮未读取或枚举sealed资料。当前批准不支持泛化、模型修复成功或生产可用结论。

实际拟合、参数优化、TSCBOP/Ag势修改、MACE训练、DFT、MD及TTM均未授权。E可依据新增用途映射整理拟合准备清单，但实际拟合仍须用户单独授权。

## 文书差异

申请MD末尾的C来源提交残留占位符 `{c}`。申请JSON提供正确来源commit ac6d98952fecc33c2de6ec952eb0ec5435234e63，相关文件SHA已核验，身份没有歧义。此为非阻塞文书问题，在本记录解释，不改写E原件。
