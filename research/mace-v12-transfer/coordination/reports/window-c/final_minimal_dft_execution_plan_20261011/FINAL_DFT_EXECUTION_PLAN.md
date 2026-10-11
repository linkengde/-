# FINAL_DFT_EXECUTION_PLAN

**计划状态：数据角色清单冻结；DFT执行仍为HOLD。** 本文件固定首轮最小数据规划，不启动DFT、MACE训练或修改势函数。

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
| C-D-AGL-1250K-rho097-seed01 | planned trajectory parent PF-AGL-1250K-rho097-seed01；source family当前UNASSIGNED | TRAIN或DEV候选，永不TEST | 等D QA；step3000结果PRELIQUID，不是液态标签 |
| C-D-AGL-1250K-rho103-seed02 | planned trajectory parent PF-AGL-1250K-rho103-seed02；source family当前UNASSIGNED | TRAIN或DEV候选，永不TEST | 真实结构及QA待核 |
| C-E-AGL-1250K-rho097-seed02 | planned trajectory parent PF-AGL-1250K-rho097-seed02；source family当前UNASSIGNED | TRAIN或DEV候选，永不TEST | 真实结构及QA待核 |
| C-E-AGL-1250K-rho100-seed01 | planned trajectory parent PF-AGL-1250K-rho100-seed01；source family当前UNASSIGNED | TRAIN或DEV候选，永不TEST | 真实结构及QA待核 |
| C-F-AGL-1250K-rho100-seed02 | planned trajectory parent PF-AGL-1250K-rho100-seed02；source family当前UNASSIGNED | TRAIN或DEV候选，永不TEST | 真实结构及QA待核 |
| C-F-AGL-1250K-rho103-seed01 | planned trajectory parent PF-AGL-1250K-rho103-seed01；source family当前UNASSIGNED | TRAIN或DEV候选，永不TEST | 真实结构及QA待核 |

通过QA后，DFT每点保存统一口径能量、64×3力；液态周期cell的应力仅在方法收敛且应力数值误差合格后用于stress训练。EAM采样的能量/力/压力不能当DFT参考。同一轨迹、restart和派生帧作为一个source family整体分到TRAIN或DEV；不能按随机帧拆分。不同seed不自动构成独立来源。TEST由独立evaluator另行持有。

参考成本：原五云端规划估计64原子液态Ag约3–8小时/DFT点，六点为18–48个四rank串行小时；六点加pilot mesh对照的阶段估计为21–64小时。均为旧规划预估，不是当前运行许可或实测保证；D结构QA通过后再按pilot实际耗时复核。

## 五、最终TRAIN / DEV / TEST划分

| Split | 固定成员 | 状态 |
|---|---|---|
| TRAIN | 原30帧保留原训练身份；本计划8点作为secondary TRAIN replay；未来液态QA和DFT通过的完整source family可加入 | 训练未授权；去重、manifest及液态QA待办 |
| DEV | 新增6个固态DFT结构；未来液态DEV须来自TRAIN之外的完整轨迹父族 | 六点来源未选、输入未生成、DFT HOLD；液态DEV未建立 |
| TEST | 仅独立evaluator保管的sealed家族；模型和阈值冻结后一次评估 | 保持封存；C不读取或枚举 |

任何同一父族的位移点、超胞、seed续段、轨迹帧、restart和数值变体整体放进同一split。公开模型评分结果只能作为回归面板，不能转成盲测。

## 停止状态

数据角色计划已冻结，但新固态来源未通过审查，E的执行准备目录未找到，所有新增DFT仍为HOLD。液态结构等待D QA。当前未启动DFT、MACE训练或修改势函数；D现有工作未被修改。
