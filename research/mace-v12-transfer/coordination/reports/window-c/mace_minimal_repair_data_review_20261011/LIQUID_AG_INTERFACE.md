# 液态Ag数据接入规范

D继续现有采样，不暂停、不修改参数。本次未发新任务。当前step3000为PRELIQUID且仍处1600K，不进入训练或DFT标签包。

候选先满足原1250K目标、热化/结构/动力学QA、安全检查及远端SHA恢复要求。候选筛选在DFT标签前完成，保留完整trajectory lineage、原子ID、组成、cell/PBC、密度/温度、seed、轨迹时间、采样器与配置SHA及选择理由。一个候选只是一个待标记几何。EAM采样能量/力/压力不得进入DFT参考字段。

后续经单独授权的GPAW/PBE单点须保存同次native extrapolated energy、force-consistent free_energy、64×3 forces及周期体相stress，eV、eV/Å、eV/Å³；stress注明ASE Voigt顺序xx yy zz yz xz xy、正负号、体积及是否电子/离子项，DFT应力不得混入采样器virial/动能压力。完整tensor与Voigt互验、应力数值收敛另审。液体拟合模型保守势对应的DFT应力，不使用MD瞬时总压力代替。

SCF、PAW/PBE、500eV、0.10eV展宽、API1e-5口径与k网格数值证据须一致；应力对cutoff/k点/Pulay效应比力更敏感，未合格应力不进入stress loss。E/F合格而stress缺失可申请energy-force候选，不声称支持R30+S。不要为补应力擅自重算已有单点。

历史R30复现保持旧native energy key；新增正式高温包建议全包统一force-consistent free_energy+同次forces，并审查所有历史帧的可追溯性。不能只对六个液态点换能量key。E0元素基准/组分和总能每原子口径要统一；不同模型零点不能直接比较。

先冻结父族和split，再取得/评分标签。同轨迹的多个帧、同初始结构的续算或冷热分支、换seed/密度后的同源后代保持同split。不同seed不自动独立。若六条任务的祖先来源审查能证明三个互不派生且非近重复家族，可提议4TRAIN+1DEV+1TEST（须按家族实际大小调整），仅为待审分配；若全部属于同一来源家族，六个点只能放一个split，不强行按帧4/1/1切。D这一条轨迹不能同时供应TRAIN和独立DEV/TEST。

C只接收TRAIN/公开DEV；TEST坐标和标签由evaluator保管，仅回传盲态家族隔离/去重证明。开发诊断过的结构不可重新标独立TEST。独立液态家族不足时可以训练候选研究，但独立液态验证保持HOLD；后续需要独立来源，不把本轮六点数目承诺为足够。

接口manifest必须有source/archive/input/result/method/PAW/config SHA、SCF状态、能量key、forces/stress units、parent/split/role审批引用及QA状态。只有角色明确授权且标签合格后才组包。本轮不复制标签、不改变角色。
