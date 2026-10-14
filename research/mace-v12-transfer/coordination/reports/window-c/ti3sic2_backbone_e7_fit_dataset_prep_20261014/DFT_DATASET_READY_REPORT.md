# Stage E7 — Ti₃SiC₂ 骨架势拟合前数据审计

## 决定

```text
STAGE_C_GEOMETRY_IMPORT=PASS
STAGE_C_DFT_LABEL_IMPORT=NO_LABELS_PRESENT
GEOMETRY_AUDIT=PASS (10/10)
FIT_READY_LABELLED_STRUCTURES=0
TRAIN_ROLE_FROZEN=0
INDEPENDENT_DEV_FAMILY=ABSENT
INDEPENDENT_TEST_FAMILY=ABSENT_NOT_ACCESSED
START_TRIAL_FIT=NO
```

本次核验纠正一个关键数据描述：Stage C 提交 `0a8585e03522b88075e5f62770a8b5409031bf69` 中的 10 项是**DFT 候选几何**，不是已经完成 DFT 的数据。Stage C 的 `FINAL_DECISION.md` 明确 `START_DFT=NO`；其十份 extxyz 只有坐标/结构元数据，不含能量、力或应力标签。因而这些几何可进入待计算清单，不能进入拟合数据集。

## 输入与结构核验

从 Stage C 分支 `cloud-c-tisc-dft-extension-plan-20261014` 导入了 `DFT_STRUCTURE_LIST.csv`、`geometry_validation.json` 和 10 个几何文件，保留原始字节并记录在 `stage_c_import/`。Stage C 原目录的 `SHA256SUMS.txt` 全部通过；逐条对比清单 SHA 后，10/10 文件哈希匹配。

ASE 3.29.0 只用于读取结构。核验结果见 `STRUCTURE_AUDIT.csv`：所有结构的原子数、Ti/Si/C 元素计数、`local_new_id` 唯一性、坐标/晶胞有限性、PBC、源端记录的最短周期距离及结构 SHA 均通过。十个最短原子间距范围为 2.074173–2.102018 Å。四个 100 原子表面结构为 TTT PBC，c 方向含 20 Å 真空；其不同终止化学计量的总能量不能直接当作可比表面能。

核验到的几何分组：

| 结构 | 数量 | 原子数/组成 | 几何来源与限制 |
|---|---:|---|---|
| 体积/剪切 | 4 | 48 / Ti24Si8C16 | 同一 Kisi/COD 基线父结构的体积 ±2%、剪切 ±0.5% 后代 |
| Si 位移 | 2 | 48 / Ti24Si8C16 | 同一父结构的 Si `local_new_id=7` 沿 z 位移 ±0.02 Å |
| clean surface | 4 | 100 / Ti52Si16C32 或 Ti48Si20C32 | 同一基线构建的 Ti/Si 终止对称 slab；参考态及顶层原子 +0.02 Å |

这 10 项全属于 `TSC_KISI_COD9009647_BASELINE_DERIVED_FAMILY`。来源不同的论文年份、终止面或扰动方式不会使它们自动成为独立父族。

## 可复用 DFT 标签清点

E5 来源提交 `195436e25d8bf8e4a524ce035b2c35c519f0bdc9` 记录了 6 个已核验的 DFT 力标签：5 个 48 原子基线/ Ti4f、C4f 位移 bulk 结构和 1 个从 Ag 界面派生的 114 原子 TSC-only slab。B 的最终结果提交 `34b8be988015af96d3fd8791fc041f6cc811e2b7` 对 bulk baseline 与四个位移点给出了已收敛能量、自由能和力；其方法清单 SHA 为 `702df71e586a61cd5785afc62653d7135a707032c91b72842b035b62198fd709`。E6 明确指出现有 bulk 记录没有可用 stress 标签；114 原子 slab 在本轮只确认 E5 力响应来源，其他标签完备性不作推断。

本次对 B 的五个 48 原子 extxyz 输出独立复核了逐文件 SHA、原子数、有限逐原子力、native energy 与结果 JSON 一致性、自由能记录和 SCF 收敛；五项全部通过。机器记录见 `B_EXISTING_LABEL_VERIFICATION.json`。114 原子 slab 的原始 DFT label 文件没有纳入 Stage E7 包，因此保留为 E5 已核验的力诊断来源，不把它计入本次独立文件校验完成的可训练数据。

这些历史标签已经用于 E5 响应矩阵与参数诊断，当前数值控制/诊断角色保持不变。它们至多是未来 TRAIN 候选，必须经独立数据用途冻结后才能回放；不得再用作 DEV 或 TEST。E5 的六个结构并不是六个独立验证样本：五个 bulk 点共用同一父结构，slab 又是界面派生结构。

| 当前库存 | 数量 | 能否立即拟合 |
|---|---:|---|
| Stage C 几何候选 | 10 | 否，缺 DFT 能量/力/应力 |
| 既有 E5 DFT 力标签 | 6 | 否，角色仍为数值控制/诊断，未冻结为 TRAIN |
| 已冻结 TRAIN | 0 | 否 |
| 独立 DEV | 0 | 否，来源家族尚未建立 |
| 独立 TEST | 0（本轮未读取/枚举封存内容） | 否，须由 evaluator 单独保管 |

`DFT_LABEL_INVENTORY.csv` 按结构列明标签组件和角色状态。`fit-ready` 数量按“有可核验标签且已获准 TRAIN 角色”计算，当前为 **0**；另有 6 个已存在的 force-labelled 诊断记录可在后续申请 TRAIN 回放。

## TRAIN / DEV / TEST 建议冻结规则

- **TRAIN 候选：** Stage C 的 10 个几何及原来 baseline 家族的全部后代只能整体留在一个组。Stage C 建议 TRAIN 候选，但没有冻结角色；完成合格 DFT 标签并取得角色批准前保持 `GEOMETRY_ONLY`。既有五个 bulk numerical-control 标签保持原角色，可另行申请整族转为 TRAIN 回放。E5 的 114 原子界面派生 slab 不得作为独立验证，可保持诊断或另行批准作为 TRAIN。
- **DEV：** 当前无可用结构。先寻找一个真正独立来源的 parent family，再为其构建参考/响应子结构；不能从现有 baseline 简单位移、扩胞、改变真空或重命名获得。
- **TEST：** 保留独立 evaluator 管理。本审计未读、未列举、未访问封存 TEST 结构或标签。TEST 父族必须与 TRAIN 和 DEV 的来源根、弛豫起点及所有后代完全隔离。

当前实际冻结数据角色计数均为零。`DATASET_ROLE_LEDGER.csv` 保存的是现有状态及候选建议，没有更改既有标签角色。

## E5 参数环境与 Stage C 候选的覆盖

E5 参数映射来源为 E5 原始势文件 SHA-256 `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`。径向外截断 `R+D` 为 Ti–C 2.918730 Å、Ti–Si 3.554301 Å、Si–C 3.211505 Å。`PARAMETER_ENVIRONMENT_ACTIVATION.csv` 给出十个 Stage C 几何逐个统计的周期 pair 数；这是几何 cutoff 覆盖，不是新的 TSCBOP 能量/力计算，也不代表参数可辨识。

- 六个 48 原子 bulk 几何各含 96 个 Ti–C 与 48 个 Ti–Si 周期 pair 进入外截断。
- 四个 100 原子 slab 几何各含 192 个 Ti–C 与 96 个 Ti–Si 周期 pair 进入外截断。
- 十个几何均没有 Si–C pair 进入 3.211505 Å 外截断；Si–C 参数 `A_Si-C/B_Si-C` 仍无支持，按 E6 建议冻结。
- 以上结构能提供 Ti–C、Ti–Si 邻域，但 pair 数不构成对 `lambda2_Ti-C` 或任何参数的力响应/秩证明。顶层位移和 shear 是否实质增加可辨识度，须等参考 DFT 标签就绪后再做独立诊断。

| 参数块 | 10 个 Stage C 几何中的当前证据 | E5 实际力响应证据 |
|---|---|---|
| `A_Ti-C` | 每个结构都有 Ti–C pair 落在外截断内 | E5 六个已标签诊断结构均有非零响应 |
| `B_Ti-Si` | 每个结构都有 Ti–Si pair 落在外截断内 | E5 六个结构均有非零响应 |
| `A_Ti-Si` | 每个结构都有 Ti–Si pair 落在外截断内 | E5 六个结构均有非零响应；与 `B_Ti-Si` 高度相关 |
| `lambda2_Ti-C` | Ti–C 邻域存在，但几何 pair 计数不能单独证明角向导数强度 | E5 六个结构均有非零有限差分响应；这仍是旧诊断，不是新增标签或拟合 |
| `A_Si-C`, `B_Si-C` | 十个结构均无 Si–C pair 进入外截断 | E5 两列严格零响应；保持冻结 |

逐结构的几何 pair 计数在 `PARAMETER_ENVIRONMENT_ACTIVATION.csv`；逐结构、逐参数的 E5 数值响应副本在 `e5_reference/force_response_by_structure.csv`。两个表所表示的量不同，不能将 cutoff 覆盖计数称作参数敏感度。

E5 实际 force-response 数据见随包只读副本 `e5_reference/force_response_by_structure.csv`，源 SHA-256 `d2b9ab6f69452c051add33708210a4f054d85b56b3ef858a0b5d04c6cdd722c8`；Si–C 截断核对副本见 `e5_reference/pair_activation_analysis.csv`，源 SHA-256 `f2030481ca2a315ef1843f7aa2366bff79ef3f0d44410fdca7839880fec74ff9`。这些是旧模型有限差分响应诊断，不能充作新 DFT 标签或拟合成功证据。

E5 全矩阵是 1062×8、rank 6/8；两列 Si–C 响应为零。剩余六列虽代数满秩但归一化条件数 1104.743，存在明显共线性。E6 给出的最小探索候选为 `A_Ti-C` 与 `B_Ti-Si` 两个 tied scalar block；只有在未来 TRAIN、独立 DEV 均准备好后，才可另行申请试拟合。`lambda2_Ti-C` 是后续三参数备选，不是本次启动许可；`A_Ti-Si` 与 `B_Ti-Si` 高度反相关，不建议与其它块一起盲目释放。Si–C 两块继续冻结。

当前请求列出的四个参数 `A_Ti-C, B_Ti-Si, lambda2_Ti-C, A_Ti-Si` 是一个拟议子集，不是 E5 所有已响应参数，也不是 E6 推荐的首轮放开组合。E5 另有 `B_Ti-C`、`lambda1_Ti-C` 两个非零响应列；E6 只建议先评估 `A_Ti-C, B_Ti-Si`，再按独立 DEV 结果决定是否加 `lambda2_Ti-C`。本报告不选择参数值或启动拟合。

## START_TRIAL_FIT 门槛

本数据包未达到启动条件。后续至少需要：

1. 先获得 10 个 Stage C 几何的授权 DFT 能量/力标签；体积/剪切点补齐一致的应力标签并确认能量口径、GPAW/PBE 设置和收敛证据。
2. 通过单独数据用途审批，将要回放的数值控制 family 整体冻结为 TRAIN；同一父族及全部派生几何不得拆分到验证集。
3. 建立且验证独立 DEV parent family 和其 DFT 标签；TEST 保持 evaluator 独立，完全不用于参数选择。
4. 冻结一个有限参数子集及训练目标/权重/约束，先从 E6 推荐二参数 `A_Ti-C, B_Ti-Si` 设计候选，不把用户列出的四块默认全放开。
5. 审核统一 DFT energy/force/stress 口径、输入/结果 SHA、原子 ID 映射和全家族来源链，并完成最终角色清单审计。

**结论：** 结构导入及几何审计已完成，但 Stage C DFT 标签不存在、DEV/TEST 独立族缺失、TRAIN 角色未冻结。拟合准备集当前为 `DATASET_READY=NO`；`START_TRIAL_FIT=NO`。没有修改 TSCBOP、运行参数优化、生成候选势或运行 MD。
