# 数据身份、父族、去重与切分规则

## 这 30 个候选的预定身份

`MASTER_TASK_MATRIX.csv` 给出唯一的提案 task ID、所属云端、建议 seed、温度/密度或响应类别、父族、atom count、来源状态和 DFT 预算。矩阵中的 parent family ID、seed 和构型模式是规划标识；结构生成/来源冻结之前不能当作已验证几何身份。全表输入状态均为 `STRUCTURE_NOT_READY`，input SHA 状态为 `NOT_COMPUTED`。

- 液态 Ag：18 个独立 trajectory seed parent，64 atoms，1050/1250/1600 K × 0.97/1.00/1.03 温度对应的液态参考密度 × 2 seeds；每个 seed 的采样轨迹与所有派生帧留在同一个 parent。1250 K 六点为这 18 行中的子集。
- Ti₃SiC₂：2 个全新来源 parent family × 每组基准及两个小热振动/模式位移，共 6 点；每个 parent 的三种配置始终同组。
- 固态 Ag：2 个全新来源 parent family × 每组基准及两种扰动/小体积响应，共 6 点；同 parent 的所有响应配置始终同组。

对于固相，表中 task-specific seed 是拟议生成/扰动 seed；source parent 仍需在实际 CIF/supercell/seed 核验后冻结。基准几何可以没有随机位移，但仍绑定源生成选择和独立 parent seed。不得根据空矩阵直接编造 atom coordinates。

## 来源与重复检查

每个未来输入需要不可变的 provenance manifest：原始 CIF/几何或 trajectory 来源、采样器和版本/许可、源文件 SHA-256、parent/seed/frame index、坐标转换、元素与固定 atom ID、cell/PBC、温度/密度或位移模式、输入字节 SHA。液态 Ag 还需采样时间、稳定取样窗、采样序列 SHA、统计量及抽帧逻辑。

冻结几何后按以下顺序去重：

1. 原始来源文件 SHA 与生成轨迹/帧 lineage；
2. 元素种类和顺序、atom ID、晶胞、PBC、坐标的精确 canonical hash；
3. 周期最小镜像下的标准化几何 hash 和最近邻距离；
4. 几何相似度/局部环境（RDF、CN/SOAP 等）近重复检查；
5. 与已公开非封存 archive、V18 train/dev、三模型已评分 parent 和 A/B 已发布参考按安全角色清单核对。

若来源缺失、SHA 冲突、与其他云端/历史数据重复或父族不清楚，标 `HOLD_DUPLICATE_OR_UNKNOWN_SOURCE`，不要静默换 ID、改 parent 或重复算。mesh、sigma、cell 数值对照是同一个几何的数值控制，不能计为独立结构。

## 五云端一致性

C 只分配唯一 task ID 和审查主矩阵。D–H owner 各用独立 branch 和本机独立的 `${DFT_WORK_ROOT}/<task_id>/`，输入/checkpoint/output/logs 不共享。需要在多台运行同一个冻结结构时，通过只读 manifest 镜像同一份输入并逐台核 SHA；每台 checkpoint 必须独立。运行后各 owner 回交结果 SHA、日志/summary、SCF/方法状态和耗时。C 发现 task/parent/input hash 重复时先冻结重复候选，不覆盖已在运行的 checkpoint。

## 数据角色

- 当前 15 个 Ag、10 个 Ti₃SiC₂ 纯相标签全部继续是 `numerical_pure_phase_control_not_training_or_independent_validation`。SCF PASS、哈希 PASS 或出现在报告里不改变角色。本轮不复制标签进训练集。
- 这 30 个新标签的 `proposed_role` 仅为“单独批准且 source/energy/force 核验之后的训练候选”；不自动具备训练资格，不能用作独立 validation/test。
- 六个 liquid pilot 的目的是核验来源、液态构型、k 点和 DFT 可行性；它们仍属于 18 个液态训练候选的预定集合，不是额外六个结构，也不是独立测试。
- 开发集/测试集必须由之后的新父族单独规划。开发和锁定测试之间、各自与训练之间均按 parent family 隔离。锁定测试的具体坐标/标签保持 evaluator 管理，C 不读不列；不将同 seed 不同帧、同晶胞位移、同界面不同 gap、mesh 对照拆到多个 split。

依据已归档 `split_plan.csv`：现有 V18 30 帧是历史训练、3 帧历史开发；测试计数不从此矩阵推导。B4 建议后续用新纯相/界面父族形成可见 development 与 evaluator-held test，但不属于此 30 点任务。任何模型、配方或参数选择前都要冻结 split manifest 和 SHA；没有独立 evaluator 时延后锁定测试，不降级成开发数据。
