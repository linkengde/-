# FINAL_DFT_EXECUTION_PLAN

**计划状态：数据角色清单冻结；DFT执行仍为HOLD。** 本文件固定首轮最小数据规划，并追加审查D窗口Ag64输出；不启动DFT、MACE训练或修改势函数。

## 来源与审查范围

C审查基线：提交 d048fdc428a1ba059e17de0cf68fcc2fa4dc432d，目录 mace_minimal_repair_data_review_20261011。E指定准备包 mace_repair_minimal_dft_preparation_20261011 在可读取的 cloud-e-b1-dft 头 39734568216e5a3d9853d7d052d28f9799172598 和当前可读 main 树中均未找到。因此，本计划的数据来源判断基于C已归档证据；E包中的输入、计算参数、资源排程及SHA尚不能独立核验。

没有改动8个旧标签的原始角色字段。下表中的TRAIN是本计划的二次使用/回放split分配，原始来源角色与审计历史必须一并保留。此冻结不启动训练。

## 一、可复用8点：分配到TRAIN回放

| structure_id | source_family / 来源 | 计划角色 | 用途 | 数据与限制 |
|---|---|---|---|---|
| Ti3SiC2_baseline_k10x10x4_v19 | Ti3SiC2_bulk_parent_v19；COD 9009647/Kisi 1998派生的48原子Ti24Si8C16 | TRAIN（replay） | 保留基础模型骨架能力；训练后仅作已知父族回归检查 | SCF/归档SHA/verifier通过；原 numerical pure phase control 角色保持不变 |
| Ti3SiC2_Ti4f_id0003_z_m020A_k10x10x4_v19 | 同一Ti3SiC2_bulk_parent_v19 | TRAIN（replay） | Ti z响应骨架回放 | 原 numerical control 角色不改；与同族所有位移一起切分 |
| Ti3SiC2_Ti4f_id0003_z_p020A_k10x10x4_v19 | 同一Ti3SiC2_bulk_parent_v19 | TRAIN（replay） | Ti z响应骨架回放 | 原 numerical control 角色不改；与同族所有位移一起切分 |
| Ti3SiC2_C4f_id0009_z_m020A_k10x10x4_v19 | 同一Ti3SiC2_bulk_parent_v19 | TRAIN（replay） | C z响应骨架回放 | 原 numerical control 角色不改；与同族所有位移一起切分 |
| Ti3SiC2_C4f_id0009_z_p020A_k10x10x4_v19 | 同一Ti3SiC2_bulk_parent_v19 | TRAIN（replay） | C z响应骨架回放 | 原 numerical control 角色不改；与同族所有位移一起切分 |
| AgSi_rigid_Ag_z_parent | A_AgSi_COD_derived_52atom_family；COD 9009647/Kisi 1998派生52原子Si终止slab | TRAIN（replay） | Ag–Si局部接触与Ag层z广义力回放 | 原公开force-response diagnostic角色保持不变 |
| AgSi_rigid_Ag_z_compress020 | 同一A_AgSi_COD_derived_52atom_family | TRAIN（replay） | −0.20 Å压缩界面响应回放 | 原diagnostic角色保持不变；不能作为独立DEV/TEST |
| AgSi_rigid_Ag_z_separate020 | 同一A_AgSi_COD_derived_52atom_family | TRAIN（replay） | +0.20 Å分离界面响应回放 | 原diagnostic角色保持不变；不能作为独立DEV/TEST |

**准入判断：允许8点在未来TRAIN回放中复用，当前仅完成数据计划分配。** 训练包引用原归档，不复制并改写原始角色。新增training-use manifest记录原角色、审批、prior-scoring披露、几何和来源SHA、父族及去重结果。与历史30帧先按精确几何及来源谱系去重。若已存在相同几何，不重复计权。

这8点都不是独立DEV或TEST。回放后它们属于TRAIN和已知回归检查。COD骨架及Ag–Si父族的位移、扩胞、重写或数值参数变体不跨split。

## 二、新增DEV：最小6个单点

新增DFT最小数量为6：4个Ti/C骨架对称位移端点和2个Ag–Si距离端点。每个结构须有收敛DFT能量与全原子力。结构家族须在任何模型评分前冻结。

| structure_id（计划ID） | 来源家族 | 角色 | 计算目的 | 状态 |
|---|---|---|---|---|
| DEV-TSC-Ti4f-z-m020 | 新Ti3SiC2 DEV父族待选；不得用已评分COD9009647骨架简单位移、扩胞或重命名 | DEV | Ti位点z向−0.02 Å；与正位移点组成恢复力斜率和逐原子力验证 | 无可执行结构；DFT未授权 |
| DEV-TSC-Ti4f-z-p020 | 与Ti负位移同一新TSC父族；来源尚未选定 | DEV | Ti位点z向+0.02 Å；计算Kz并检查两端力偏置 | 无可执行结构；DFT未授权 |
| DEV-TSC-C4f-z-m020 | 与Ti位移点同一新TSC父族，除非来源登记另行定义 | DEV | C位点z向−0.02 Å；测试C骨架恢复力 | 无可执行结构；DFT未授权 |
| DEV-TSC-C4f-z-p020 | 与Ti位移点同一新TSC父族，除非来源登记另行定义 | DEV | C位点z向+0.02 Å；计算C位点Kz及力偏置 | 无可执行结构；DFT未授权 |
| DEV-AGSI-gap-m010 | 新Ag–Si界面父族待选；旧COD/A的52原子家族不合格作独立源 | DEV | 固定合法registry下d0−0.10 Å，测界面接触力 | 来源、termination、registry未冻结；DFT未授权 |
| DEV-AGSI-gap-p010 | 与gap负端点同一新Ag–Si界面父族 | DEV | 固定同一registry下d0+0.10 Å，测Ag层广义力和局部恢复响应 | 来源、termination、registry未冻结；DFT未授权 |

这6点是最低局部诊断屏，不是统计充分性的保证。四个TSC位移共享一个父族；两个Ag–Si间距共享一个父族。若没有经来源核验的新父族，就不能声称是独立DEV。三个历史TSC论文描述同一晶体原型；仅换论文、CIF编号或晶格参数并不能证明独立。Zhang 2026的黏附表没有提供可复现坐标和逐原子力。结构须在DFT前完成来源、完整化学组成、cell/PBC、原子映射和周期近重复审查。

## 三、最低DFT方法和资源预算

单点输出能量和完整原子力。建议PBE/PAW、PW cutoff 500 eV、展宽0.10 eV、SCF阈值1e-5；所有设置仍须随来源结构核准。TSC候选k点8×8×4与10×10×4同几何对照；slab候选6×6×1与8×8×1、xy偶极修正，真空至少18–20 Å并对同几何增加4 Å检查。不能未经测试一律设Gamma。它们是候选收敛方案，不是已经收敛的计算输入。

可在4个周期TSC点中同时输出bulk应力，供stress路线开发；只有应力通过cutoff和k点收敛审查才可用。含真空界面的3D应力不用于应力监督。

资源估算按4 MPI rank：4个48原子TSC单点约3–6小时/点，12–24小时；2个80–120原子slab单点约6–12小时/点，12–24小时；总计约24–48四rank串行墙钟小时（96–192 rank-hours）。不含k网格/真空对照、排队、失败重跑和结构来源准备，实际可能更长。E的具体执行预算尚无法审阅，因为其准备目录当前不可见。

## 四、液态Ag接口：等待D QA

D/E/F的六条Ag64任务只可成为TRAIN或DEV候选，不能进入TEST。D step3000当前分类为PRELIQUID，仍处1600 K热段，不具备1250 K液态QA资格，也尚无可用DFT标签。

| structure_id | source_family | 角色 | 用途与当前状态 |
|---|---|---|---|
| C-D-AGL-1250K-rho097-seed01 | trajectory parent/source family `PF-AGL-1250K-rho097-seed01`；共同source root为六任务campaign | **HOLD；step3000和step8000同族，当前不进入TRAIN/DEV/TEST** | step3000为PRELIQUID；step8000 QA数值通过但包内`STRUCTURE_PASS=F`冲突，待纠正后再考虑DEV |
| C-D-AGL-1250K-rho103-seed02 | planned trajectory parent PF-AGL-1250K-rho103-seed02；source family当前UNASSIGNED | TRAIN或DEV候选，永不TEST | 真实结构及QA待核 |
| C-E-AGL-1250K-rho097-seed02 | planned trajectory parent PF-AGL-1250K-rho097-seed02；source family当前UNASSIGNED | TRAIN或DEV候选，永不TEST | 真实结构及QA待核 |
| C-E-AGL-1250K-rho100-seed01 | planned trajectory parent PF-AGL-1250K-rho100-seed01；source family当前UNASSIGNED | TRAIN或DEV候选，永不TEST | 真实结构及QA待核 |
| C-F-AGL-1250K-rho100-seed02 | planned trajectory parent PF-AGL-1250K-rho100-seed02；source family当前UNASSIGNED | TRAIN或DEV候选，永不TEST | 真实结构及QA待核 |
| C-F-AGL-1250K-rho103-seed01 | planned trajectory parent PF-AGL-1250K-rho103-seed01；source family当前UNASSIGNED | TRAIN或DEV候选，永不TEST | 真实结构及QA待核 |

### D窗口新增Ag64结构的角色复核（2026-10-11）

本次核对D分支 `cloud-d-b1-dft` 的最新可见提交 `2ad687d7f74e035a9d2274e7d1e9a70577eaa76f`。该提交只包含任务 `C-D-AGL-1250K-rho097-seed01` 至step3000的分段归档；没有发现已冻结的液态候选结构文件或其DFT标签。step3000归档SHA-256为 `cd00e213b89242fa9501aacbadd5cec8c23e938a61aeebbb2eb138997deea7d2`，与远端sidecar及receipt一致。C已有离线审查记录将其判为 `PRELIQUID`：结构和移动性呈热态无序迹象，但仍处1600 K阶段，温度末窗未稳定，也未到1250 K；source manifest的 `STRUCTURE_PASS=false`、`liquid_QA=NOT_PERFORMED`、`DFT_LABELS=ABSENT`。因此它目前是诊断帧，不是经QA确认的液态Ag数据点。

| structure_id | source_family | role | reason |
|---|---|---|---|
| `C-D-AGL-1250K-rho097-seed01-step3000` | 计划trajectory parent为 `PF-AGL-1250K-rho097-seed01`；实际source family仍未确认；祖先source root为 `PF-WMH2006-Ag64-NVT-melt-cool-1250K-six-pilot-source-campaign-20261011` | **HOLD_PRELIQUID（不属于TRAIN/DEV/TEST split）** | 1600 K热段且热化未建立；液态QA未通过；没有DFT能量/力标签。即使以后证明同seed轨迹连续，也必须把该轨迹及全部restart/派生帧作为一个parent family，不可跨split。 |

**TRAIN/DEV决策规则：** 液态结构须先完成液态QA并取得合格、统一口径的DFT能量与全原子力（应力仅在数值收敛后使用）。此后按完整 `source_family` 分配：若该来源族被纳入拟合或纯相回放，则整个族进入TRAIN；若为评估液态泛化而预先留出的独立来源族，则整个族进入DEV，仅用于预先定义的验证与模型选择，不进入拟合。source family未确认、与TRAIN存在共同结构/轨迹祖先、或已用于模型选择评分的家族不得称为独立DEV。不得按同一轨迹中的帧、不同续算段或随机seed切分一个家族。D当前step3000帧不满足TRAIN或DEV准入条件，不能仅因它来自D窗口而默认入组。

**TEST隔离：** 以上D数据不得进入TEST，也不改变任何sealed TEST的归属或访问规则。液态独立TEST继续由独立evaluator管理；C不读取其结构或标签。

### Ag64 step8000候选交接的角色审核（D提交38ff5afb）

审核对象为D分支提交 `38ff5afb5c3cbe14ca0889212910f0ca6558e991` 的 `Ag64_step8000_candidate_release/`。`SHA256SUMS.txt`中四项哈希均与Git远端文件一致；MANIFEST中结构、QA和热日志哈希也匹配。结构为64个唯一原子ID的Ag，立方晶胞边长10.826824658074 Å，三向PBC。冻结势函数SHA为 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`，与项目指定值相同。

来源核对：task ID为 `C-D-AGL-1250K-rho097-seed01`；原始source manifest给出trajectory parent/source family `PF-AGL-1250K-rho097-seed01`、campaign root `PF-WMH2006-Ag64-NVT-melt-cool-1250K-six-pilot-source-campaign-20261011`和PCG64 seed `125009701`。step3000的PRELIQUID帧与step8000属于同一轨迹族，不能拆分到不同split。D候选包没有`source_family`字段；本计划按上游source manifest中的`trajectory_parent`登记家族。不同任务seed本身不能证明D/E/F家族彼此独立，其他家族仍须来源和几何去重审查。MANIFEST/QA的step8000 checkpoint及state/array哈希彼此一致，但checkpoint本体不在交接包或当前D分支树内，因此C核实的是哈希记录一致性，尚未对checkpoint实际字节重新取回验哈希。

QA的数值结论支持“1600 K热态液体候选”：step5000和6000热均温门槛未过；step7000、8000两个连续窗口均报告三项标志全通过。两窗口均温为1637.23 K、1638.76 K，标准差155.98 K、84.68 K（目标1600 K、均温容差50 K、标准差门槛329.18 K）。局域Q6均值为0.37675、0.39634；step8000窗口的Q6斜率为+0.04510/1000步，但双侧p=0.08264，按该QA预设规则通过。配位数均值为10.5104、9.9896；每个窗口内半窗漂移0.125、0.104，均通过，但相邻窗口均值下降约0.521，显示结构统计仍有变化。MSD斜率为0.002201、0.001801 Å²/fs，拟合R²为0.9919、0.9887，均为正且通过；晶体样连通簇比例为0，最小间距分别为2.207、2.314 Å。上述筛查基于有限的250 fs窗口；QA的液态结构指标使用局域Q6加相干键图，非CNA/PTM。温度目标为1600 K，所以这不是已平衡的1250 K样本。

**当前结论：`STRUCTURE_STATUS=HOLD`。** MANIFEST和QA把结果记为`TWO_CONSECUTIVE_WINDOWS_ALL_VALID`，而冻结XYZ注释同时写有`candidate_status=PRELIMINARY_LIQUID_CANDIDATE`及`STRUCTURE_PASS=F`。在字段含义未由D澄清、或D未发布保留原件的更正/补充记录前，不能把false标记视作无害旧元数据。加上checkpoint本体当前不可独立取回验证，本轮不批准这帧进入DFT。

**条件性数据角色建议：DEV候选，不建议TRAIN。** 解决标记冲突、核验step8000 checkpoint实际内容并完成与其他液态轨迹的source-family/几何去重后，可把结构 `Ag64_step8000` 作为固定来源族 `PF-AGL-1250K-rho097-seed01` 的独立液态开发验证候选。理由是当前没有可用的液态Ag DFT训练标签，保留该完整来源族用于比较模型在液态环境的误差更有验证价值；一旦用作DEV评分，不得再转为TRAIN。液态TRAIN数据应来自另一条独立、去重通过的完整来源族。该建议不等于DFT启动许可。此族不进入TEST。

| structure_id | 当前role | 条件性DFT role | source_family | DFT类型与标签需求 |
|---|---|---|---|---|
| `Ag64_step8000` | **HOLD** | **DEV**（仅在冲突和来源核验完成后） | `PF-AGL-1250K-rho097-seed01`；source root为六任务campaign | 固定几何、固定晶胞的PBE/PAW周期性单点；要求统一口径的native energy、force-consistent free energy及64×3全原子力。输出完整应力张量供审查；只有k点/应力收敛审查通过后才可用作stress监督。采样EAM的能量、力、压力不是DFT标签。 |

DFT标签角色一旦批准为DEV，全家族的可选帧和同源续算都留在DEV；step3000的PRELIQUID帧不作为标签，但同样不得迁入TRAIN。sealed TEST继续由独立evaluator保管，C不读取、不枚举。

通过QA后，DFT每点保存统一口径能量、64×3力；液态周期cell的应力仅在方法收敛且应力数值误差合格后用于stress训练。EAM采样的能量/力/压力不能当DFT参考。同一轨迹、restart和派生帧作为一个source family整体分到TRAIN或DEV；不能按随机帧拆分。不同seed不自动构成独立来源。TEST由独立evaluator另行持有。

参考成本：原五云端规划估计64原子液态Ag约3–8小时/DFT点，六点为18–48个四rank串行小时；六点加pilot mesh对照的阶段估计为21–64小时。均为旧规划预估，不是当前运行许可或实测保证；D结构QA通过后再按pilot实际耗时复核。

## 五、最终TRAIN / DEV / TEST划分

| Split | 固定成员 | 状态 |
|---|---|---|
| TRAIN | 原30帧保留原训练身份；本计划8点作为secondary TRAIN replay；未来通过QA及DFT标签核验、且被选作拟合的完整液态source family可申请加入 | 训练未授权；没有可用液态TRAIN标签；需要其他独立轨迹族 |
| DEV | 新增6个固态DFT结构；`PF-AGL-1250K-rho097-seed01`为条件性液态DEV候选 | D Ag64 step8000当前HOLD；纠正结构标志、验证checkpoint与source-family去重前不进入DEV/DFT |
| TEST | 仅独立evaluator保管的sealed家族；模型和阈值冻结后一次评估 | 保持封存；C不读取或枚举 |

任何同一父族的位移点、超胞、seed续段、轨迹帧、restart和数值变体整体放进同一split。公开模型评分结果只能作为回归面板，不能转成盲测。

## 停止状态

数据角色计划已更新，但新固态来源未通过审查，E的执行准备目录未找到，所有新增DFT仍为HOLD。D的step8000 QA窗口数值通过，然而XYZ的`STRUCTURE_PASS=F`与QA/MANIFEST不一致，且checkpoint本体未能从交接包/当前D树取回，因此当前结构角色为HOLD；条件性建议整族作DEV、不作TRAIN或TEST。当前未启动DFT、MACE训练或修改势函数；D现有工作和归档未被修改。
