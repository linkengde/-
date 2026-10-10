# 云端 C：三模型冻结势静态对比

状态：**两模型已完成相同公开开发参考上的冻结推理；MATPES-PBE-0 因 ASL 许可范围未确认且官方模型元数据请求返回 HTTP 403，未下载、未推理。** 本报告使用 MACE-MP-0b3-medium 与冻结 V18 的结果，不将缺席的 MATPES 数值化或推断。所有公开参考只作开发诊断/数值控制，不是独立验证集。

## 核心结论

三模型比较尚不能完整给出数值排名，但已有对照明确显示两种已测模型各有失效域：

- 原始 MACE-MP-0b3-medium 保留了 Ti₃SiC₂ 骨架力；10×10×4 基准全原子力向量 RMSE 为 0.038173 eV/Å，而 V18 为 2.888854 eV/Å（75.7 倍）。同骨架 Ti 位移 ±0.02 Å 点的 RMSE 为 0.039–0.041 对 2.881–2.897 eV/Å；C 负位移点为 0.038849 对 2.886734 eV/Å。V18 的 Ti/C 误差主要沿 z。该结果强化了“从基础模型重新初始化”的选择，但不能证明目前基础模型已经达到生产要求。
- 六个 Ag–C/Ag–Si/Ag–Ti 局部界面开发构型上，V18 的全原子向量 RMSE 低于基础模型：按接触分别为 0.141/0.057/0.104 与 0.906/0.736/0.519 eV/Å。两者都未达到暂定 0.05 eV/Å 门槛。局部标记原子对分离力绝对误差均值却是基础模型较低：Ag–C 0.108 对 0.217，Ag–Si 0.098 对 0.188，Ag–Ti 0.200 对 0.232 eV/Å。整体力向量分数不能代替目标接触力检查。这六点源自已复用的开发父族，并非独立验证。
- 纯 Ag 七个基准/小位移/剪切/体积结构的 pooled force-vector RMSE 为基础模型 0.000648、V18 0.001271 eV/Å。一个完美晶格近零力点只说明该结构的静态对称性；不等于缺陷、热振动或熔化能力通过。Ag 单原子对称位移恢复斜率（DFT/基础/V18）为 5.863/5.302/4.747 eV/Å²。三点 Ag 剪切能量曲率为 174.0/174.2/164.9 eV/cell；log-volume 曲率为 375.8/381.5/316.9 eV/cell。它们只是当前 PBE 父族和指定应变括区的有限响应回归。
- 已完成并核验的 A 端 Ag 层 −0.20 Å slab 端点上，基础模型/V18 全原子 RMSE 为 0.202/2.635 eV/Å；DFT/基础/V18 的整层 z 广义力为 −1.357/−3.381/−4.162 eV/Å。这个单点不能给出成对斜率或完整分离曲线；A 的 +0.20 Å 尚无通过核验的归档。仓库 main 的 A 任务卡仍将 −0.20 Å 显示为 running iteration 50，但该标签自带 `progress.json` 为 complete/SCF 59、summary 为收敛且归档 verifier 为 PASS；此处只依据自洽归档评分，并将 A 任务卡状态滞后记为交接记录，不触碰 A 队列。现有 Ag–Si slab 家族与六个局部界面不是相互独立的统计测试。
- 已发布的旧 EAM+Tersoff+Morse 数据另存在 `old_hybrid_vs_dft.csv`，明确注明来源为历史 LAMMPS run0。本轮没有运行 LAMMPS；不将混合势误差并入任何 MACE 汇总。

因此，下一轮获授权的 V19 受控实验优先以精确哈希的 MACE-MP-0b3-medium **重新初始化**，并以冻结 V18 作非候选基线，按既有三臂方案比较“原配方重启 / 冻结低层 / 骨架回放”。但当前纯相标签均保持 `numerical_pure_phase_control_not_training_or_independent_validation`，不得仅凭此审计把它们转为训练标签；R30+S 数据臂仍须先解决独立数据角色审查。V18 在局部界面向量误差较低不构成继续从 V18 微调的理由。

MATPES-PBE-0 值得在许可范围获确认后作为候选起点做同结构推理：官方目录注明 MATPES-PBE 与 PBE 且无 +U，理论上比 MACE-MP-0b3（目录标为 MPTrj、PBE+U）更贴近当前 PBE DFT 设置。但许可和权重阻碍未解除，因此不能声称它更好，也不把它加入 V19 执行计划的已测臂。MACE-MP 的训练理论与本仓库 PBE 标签不同，也应在后续归因中保留为方法差异。

## 数据和口径

审计了 22 个公开、未封存结构：8 个固态 Ag 结构（6×6×6/8×8×8 同 parent 基准、±0.01 Å 单原子位移、±0.5% shear、−2%/+2% volume）；5 个 Ti₃SiC₂ 结构（10×10×2/10×10×4 基准、Ti ±0.02 Å、C −0.02 Å）；6 个 Ag–X 局部接触结构（Ag–C、Ag–Si、Ag–Ti 各正负局部位移）；以及 3 个 Ag–Si slab 参考（原真空、26 Å 真空数值对照、A 已核验 −0.20 Å 刚性 Ag 层端点）。逐结构来源/结果/摘要哈希、SCF 与 verifier 状态、数据角色、父族和训练资格见 `reference_inventory.csv`。

输入仅在 `summary.scf_converged=true`、归档 verifier 全 PASS、结果/源哈希一致后纳入。基础模型与 V18 均在 CPU float64 推理，元素次序、原子 ID 顺序、PBC/cell/positions 与参考力逐项核对；八个已发布审计标签的 MACE/V18 推理逐原子力复现检查通过。力向量 RMSE 定义为 `sqrt(mean_i(||F_model,i − F_DFT,i||²))`；另报告元素、x/y/z RMSE 和最大逐原子向量误差。`per_element_errors.csv`、`per_atom_errors.csv` 保留结构与模型维度。

能量只在同一组成、同一父族、同一模型内部相对 parent 比较，不跨模型比较绝对能量零点。剪切与 log-volume 使用 parent/minus/plus 三点拟合并同时列出 parent 处一阶、二阶导数；volume 实际采用 −2%/+2% 体积括区，在 log-volume 上不是严格等距。DFT native energy 与 force-consistent free energy 分开列出。尚未解读固液共存、液态 Ag 或冷却路径物性。

已完成 Ti ID3 的 ±0.02 Å 配对恢复响应，使用零基数组索引 2、跨度 0.04 Å；恢复斜率 −ΔFz/Δz（DFT/基础/V18）为 17.471/13.868/6.292 eV/Å²。B 的 C ID9 +0.02 Å 尚在运行（本次快照来自 `main` commit `1ea61a2` 的 iteration 24 记录），故只保留 C− 的点误差，不计算 C 成对斜率。A 的 Ag layer +0.20 Å 也仍在运行，故只报告已完成的负端点力。

## 熔池与凝固适用范围

这批静态固相及局部界面参考没有液态 Ag、过冷 Ag、固液共存界面、高温 Ti₃SiC₂、热态 Ag/骨架界面或短程高温排斥标签。也没有熔点、扩散、密度、潜热、RDF、结晶序参量或界面迁移的势函数验证。因此本结果**不证明熔化、凝固、喷溅或熔坑适用性**。相关数据分组、DFT 参考和物性验证仍按已归档 V19 设计方案，须另行明确授权后执行。

## 复现与文件

冻结推理复用现有 MACE 环境（Python 3.12.14、PyTorch 2.5.1+cpu、mace-torch 0.3.16、ASE 3.29.0、NumPy 1.26.4、SciPy 1.17.1）。推理只使用四 CPU 线程，脚本见 `compare_static.py`。MATPES 未下载或加载。

查看：`MODEL_SOURCES.md`、`V19_REPAIR_PLAN_REVISION_20261010.md`、`model_inventory.csv`、`reference_inventory.csv`、`structure_metrics.csv`、`per_element_errors.csv`、`per_atom_errors.csv`、`ag_response.csv`、`ag_elastic_response.csv`、`backbone_interface_responses.csv`、`interface_separation_and_layer_forces.csv`、`old_hybrid_vs_dft.csv`、`reproduction_checks.csv`。`SHA256SUMS.txt` 给出本目录交付文件摘要。

禁止工作均未启动：DFT、训练、MD、LAMMPS、TTM；未读封存标签、改权重、改原始标签或 A/B 任务，也未将数值控制标签升级为训练用途。
