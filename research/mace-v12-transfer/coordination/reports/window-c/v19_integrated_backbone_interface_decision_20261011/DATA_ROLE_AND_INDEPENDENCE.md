# 本次核验更新

本次白名单是旧25项加已公开C+一项，共26项：15 Ag记录对应7个输入SHA，11 Ti₃SiC₂记录对应5个输入SHA；这是几何身份计数，不是独立父族计数，也不是全仓库数据总数。26项SCF、归档verifier、SHA均通过，native/free_energy口径均有记录。全部训练资格仍为NO，角色变更数为0。逐条证据见`pure_phase_role_audit.csv`和`role_audit_summary.json`。以下保留旧快照与新候选设计的区别。

# 固相修复的数据角色与独立父族审查（离线设计）

本轮仅只读审查公开 C 清单和方案，并提出待审设计；没有新增任务登记、生成结构、运行 DFT/训练/MD，也没有更改任何数据角色。D/E/F 已授权六液态 pilot 的所有权和既定流程保持；G/H 待命、其余24项HOLD保持。本文件中的 proposal_group 是审查建议分组，不是可执行 task ID。

## 读取的安全来源与证据限度

已读取 `v19_controlled_repair_design_20261010/{PLAN.md,pure_phase_source_eligibility.csv}`、`three_model_static_comparison_20261010/{reference_inventory.csv,V19_REPAIR_PLAN_REVISION_20261010.md}`、`unified_potential_independent_review_20261010/{reusable_dft_inventory.csv,V19_READINESS.md}` 及 `unified_potential_dft_data_plan_20261010/{existing_data_inventory.csv,proposed_batches.csv,split_plan.csv,METHOD_AND_ROLE_RULES.md}`。优先使用这些 C 元数据，没有打开 V18 训练30帧/历史dev3帧的坐标或标签；没有枚举/读取任何既有封存结构或标签。

旧资格快照含15 Ag与10 Ti3SiC2记录，分别只有7与4个 source_geometry/source SHA；这不是7或4个独立父族。全部 `train_eligible_now=False`，角色为 `numerical_pure_phase_control_not_training_or_independent_validation`。同一Ag父源的位移/shear/volume以及TSC父源的Ti/C位移共享家族；k网格变化更不产生新几何。

旧25记录是快照计数，不能声称就是当前归档总数。旧TSC资格清单包含baseline与Ti±/C−三端点，未列后来完成的C+；根 C 本轮已结合最终四点审计并追加 C+ 公共元数据白名单核验，完整五构型归档通过；新增角色授权仍未获得。

`reference_inventory.csv` 中22公开评分结构、旧六个Ag-X残差、V18历史dev3帧以及A/B响应父族均已被看到/评分，不能作为新的独立测试。V18历史训练30帧都含Ag，没有纯Ag或纯Ti3SiC2；原训练角色并不免除再次建集时的来源、SCF、原子身份和能量口径审查。

旧CSV有两个摘要格式缺陷：`existing_data_inventory.csv` 一个公开角色摘要行包含未引用的逗号，`proposed_batches.csv` 的B1-LIQ末尾包含未引用逗号，Python DictReader产生多余字段。没有改写这些冻结报告；不能将错位的短文本当成训练资格。有效资格以逐条pure/reference账本与METHOD/REPORT为准。本轮新CSV应严格校验每行字段数。

## 训练修复候选与验证必须分别处理

### 最低新增成本的训练候选：最多免5个同几何单点，当前可组装为0

优先提出未来对已核验完整TSC骨架baseline+Ti±0.02Å+C±0.02Å五构型的逐条角色复核。若用户以后明确批准用途变化且所有五点来源/SCF/SHA/能量口径通过，可选每个几何一个已收敛网格标签进入单独骨架回放候选包，最多避免5个相同几何的重复DFT，粗估15–30串行四rank小时（60–120rank-hours）。旧快照直接支持4个几何；根 C 本轮对白名单追加的 C+ 归档及完整五几何元数据复核已通过，五构型可供未来用途复核，但训练资格仍为0。

**当前没有角色授权，故训练可组装标签数为0，实际节省为0。** 不复制控制标签进训练目录，不在CSV把控制角色改成training。若角色将来获准转为训练，其同父源全部后代从fresh dev/test资格排除；原数值控制来源及变更记录仍永久保留。响应包覆盖的是一个父源的局部Ti/C z刚度，不能证明完整骨架构型空间足够。

不准复用时的替代提案是1个待选完整Ti3SiC2训练来源根×5构型：baseline、Ti局域z±0.02Å、C局域z±0.02Å，48–96原子，5个新增单点约15–30四rank串行小时。它也未获得结构或DFT授权；不作为当前30任务自动增项。该训练根与下述dev/test根必须分开。未来可先承诺来源/角色，再决定是否申请它，不重复既有同几何DFT。

Ag-X历史30训练帧、旧残差/界面档案只可按各自历史角色和来源审查决定是否复用；残差“training acquisition”标签不等于本轮已取得训练许可。避免把6个残差或48界面档案重算成新增标签数量。

## 最小新独立固相筛查设计：22个候选验证单点

附CSV列出待审设计，不包含真实CIF、原子坐标、已有封存ID或可运行输入。

1. **Ti3SiC2：DEV与独立TEST各1个待选来源根，每根5点，共10点。** 每根为48–96原子的完整周期Ti3SiC2晶体，组成分别可为Ti24Si8C16（48原子）或Ti48Si16C32（96原子）；baseline及Ti/C两个独立局域z位移的±0.02Å均保留在该split内。新dev源需真实独立获取/测量/构造依据，且与历史父源不存在派生关系或等价/近重复；TEST源由独立evaluator选择与保管。不能用另一个CIF编号、seed、超胞旋转、应变或复制就直接认定新根。如果找不到有证据的新根，该设计保持 `BLOCKED_SOURCE_NOT_SELECTED`。
2. **Ag–Ti、Ag–Si、Ag–C：每类DEV与TEST各1个待选真实界面来源根，每根2个固定gap点，共12点。** 每点拟80–120原子，保留完整Ti3SiC2骨架与Ag接触环境，不用孤立二原子模型替代；每根按被冻结合法界面锚点 `d0±0.10Å` 选两点（实际安全范围需源审查，未生成几何），所有gap后代同组。Ag–Ti可提出历史未覆盖的Ti暴露棱柱/台阶接触；Ag–Si可提出不同真实Si终止或重构接触；Ag–C可提出可解释的C暴露混合面/台阶局域环境。这里只是机制筛选标准，不声称表面稳定或任何这些根已独立。TEST具体来源和坐标不发送C。
3. **分配统计、物理覆盖与统计限制：** 22点分为11DEV与11未来evaluator TEST，共8个待选来源根（每split各1TSC+3接触根）。1根/split/机制仅是低成本筛查下限，不达到旧方案推荐的每机制至少2独立父族的确认强度。纯Ag新的独立dev/test及液态、凝固验证仍欠缺，本最小固相范围不新增它们、不宣称全四元素或相变生产势通过。

成本沿用公开旧规划的未实测估计：TSC3–6 h/48–96atom单点；界面6–12 h/80–120atom单点；均按4MPI ranks单作业、串行墙钟统计。10TSC为30–60 h，12界面为72–144 h，合计102–204 h（408–816rank-hours）。若另申请5新训练点，增加15–30 h，总117–234 h（468–936rank-hours）。这些不是性能保证，也不能说明22点已证明足够。

历史48atom TSC约3–4.2 h、52atom slab约3.13/3.81/4.38 h是量级依据，不能直接按原子数线性放大为80–120atom保证。大端96atom TSC或120atom slab可能超过上述6/12 h区间，尤其更密k网格、bands/SCF或真空扩尺寸改变内存/时间；原子数本轮未冻结，所有新固相运行预算均未授权，不能把液态pilot的12 h上限自动套作新固相预算。本筛查应优先在保留完整化学/表面层数且有限尺寸合格的最小真实盒审定原子数，再以获准首例重估，不以压小盒或不经验证Gamma强行满足时间预测。

同几何k点/真空等数值控制必须单列：若未来最低每个新来源根需1个更密网格对照，则2个TSC根与6个界面根可额外需要8个数值运行、42–84四rank串行小时；标签+这一示例控制为144–288 h（576–1152rank-hours）。若真空或第二个k点仍未收敛需另审预算，不自动扩算。数值控制不多算成训练构型、不产生新父族；通过的同方法同几何控制可避免重复同一个最终单点，但角色仍须另审。

若用户坚持沿用旧方案的至少2来源根/split/机制，则上述验证包应扩为20TSC+24界面=44候选验证单点，约204–408四rank串行小时（816–1632rank-hours），同几何数值控制另计。这是后续确认方案，不是本轮任务增加或已有24HOLD的自动解锁。旧B4规划48个验证点包含另外12纯Ag，TSC每根3帧；本22/44提案改为TSC每根5帧以分别覆盖Ti和C对称恢复斜率，不能把新数字写成旧批次已经获得批准。

## 来源与独立性门：先证明，再计为新父族

1. 冻结公开训练/历史评分家族账本：C旧30训练、3dev、22诊断、A/B响应、旧残差、数值控制源及全部已知祖先/转换。具体source SHA与父族关系需从允许读取的公开manifest逐例补足，不能按记录数推父族数。
2. 新DEV的source论文/样品/测量记录、源文件SHA、原子映射/拓扑、cell/PBC、termination/registry及所有生成步骤必须完整；使用独立来源不自动免除几何等价检验。不同数据库ID可能重复收录同实验，同prototype的CIF亦可能是同几何。
3. 统一origin lineage DAG。old root的k网格、真空、位移、不同seed/热帧、重新写CIF、supercell、剪切/应变和registry曲线后代保持同源分组，不随机跨split。新的真实界面机制只有在来源审查和几何/局部环境检查证明不属于旧派生/重复家族后才可接受；同一底层晶体源的共用祖先不能隐藏。
4. 精确几何SHA用于文件身份；canonical periodic symmetry/atom mapping去重及局部environment descriptor用于等价/近重复筛查。SHA不同不证明独立，descriptor阈值须在看标签前冻结。构型父族与局部环境都有差异仍不证明统计充分性，必须如实列来源层级。
5. C可验证公开训练/DEV之间的冲突，但不能查看封存TEST。独立evaluator对其持有的既有封存池和新TEST候选做盲态冲突审查，仅向C提供签署的 `DISJOINT_FROM_TRAIN_AND_DEV=YES/NO`、`DISJOINT_FROM_EXISTING_SEALED=YES/NO`、`SOURCE_LINEAGE_REVIEW=YES/NO` 及审查版本/承诺哈希；不透露封存ID、坐标、标签或冲突详情。C只能说“待独立evaluator attestation”或“收到盲态YES”，不能声称自己已读取/验证排除sealed冲突。NO或缺回执即HOLD。
6. 新TEST由evaluator预先封存，C不接收坐标/标签；模型、超参数、输出checkpoint和门槛冻结后一次评分。失败后的新一轮选择必须形成新的开发过程/测试计划，不能反复用同封存测试调参。

## DFT口径与验收

GPAW26.7、ASE3.29、PAW1.2.1、PBE PW500eV、sigma0.10eV、既有API energy/density/eigenstates各1e-5口径；max160只是本新固相设计待批准上限，不自动沿用液态pilot的授权。TSC体相PBC xyz；slab界面PBC xy/z false、审查20–26Å或更大的实际真空与z向dipole设置；k点按新cell的倒空间密度审查，禁止通用Gamma。更密网格/真空需要逐来源固定几何数值证据。

每个标签保存native extrapolated energy、force-consistent free_energy和全原子3Dforces，单位eV/eVÅ−1；固定顺序/atom IDs/cell/PBC和Ti/C位移或Ag层mask；SCF收敛、有限结果、source/input/output/PAW/method SHA及归档verifier全部PASS。有限展宽正式热训练包采用同一次free_energy+forces，不混native与free_energy；历史R30复现实验全臂同一native key，不能只换新增点口径。模型之间不同绝对能量零点不直接比较。

骨架每根baseline及Ti/C对称斜率分别报告，`Kz=-(Fz(+0.02)-Fz(-0.02))/0.04`；全原子向量RMSE、Ti/C分元素xyz和最大原子误差独立列出。界面每类报告全原子、元素xyz、接触投影、整Ag层广义力及两个gap端点中央响应；不能用总平均改善掩盖骨架或某接触类退化。待审门槛只能在DEV验证覆盖后预注册，不打开TEST调门槛。

## 当前状态

所有新根均 `UNSELECTED`，所有新单点 `NOT_AUTHORIZED/NOT_GENERATED`。旧纯相控制训练角色仍未授权，训练包不组装。本文件不登记新task、不分配到D/E/F/G/H、不解锁任何HOLD项、不读sealed；只把最低固相数据资格与来源独立性缺口明确给下一步审定。
