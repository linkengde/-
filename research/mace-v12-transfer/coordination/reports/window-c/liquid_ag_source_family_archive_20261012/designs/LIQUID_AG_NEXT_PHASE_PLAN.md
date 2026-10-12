# 液态 Ag 下一阶段数据与 DFT 申请方案

日期：2026-10-12  
负责窗口：D  仅限液态 Ag 数据准备  
状态：审查计划；未启动 DFT、MD、AIMD、LAMMPS 或 MACE 训练。

## 范围与保护项

本报告只审查 D 当前可读取的液态 Ag 轨迹、checkpoint、候选交接包和既有液态 Ag 规划记录。没有写入或重算原子结构，没有更改势函数、数据角色、Ag64 step8000 文件、E 的 Ag–Si registry DFT，也没有读取或改动封存 TEST 标签。

## 1. 已有候选和 pilot 盘点

| structure_id | source_family | 温度历史 | QA 状态 | 满足液态筛选条件 | 当前角色 |
|---|---|---|---|---|---|
| Ag64_step8000（结构 SHA-256 `08d5e9a7500e264cce2382c07d76afe032e257d5ae7f7d986868e2c646316b85`） | `C-D-AGL-1250K-rho097-seed01`；64 Ag，同一连续轨迹 | 热浴目标日志从 step 0 到 8000 均为 1600 K；记录到的实际温度范围约 1144.9–4266.3 K，早期有明显瞬态。trajectory 覆盖至 4000 fs。任务名中的“1250K”与日志目标 1600 K 不一致。 | step 7000、8000 两个连续 1000 步窗的 Thermal、Liquid_structure、MSD 标志均通过；step8000 平均温度约 1638.8 K、标准差约 84.7 K，MSD 斜率为正且拟合 R² 约 0.989。结构筛查用 Q6 与 coherent-bond 图替代 CNA/PTM。 | **部分通过**：通过本次在线液态候选筛查；不能作为已满足完整液态 Ag 验收或 1250 K 条件的证据。QA 有窗口 RDF 数组，但没有独立的 RDF 尾部/S(k) 通过判据，也没有压力/能量平稳性验收与来源族复核。 | `ARCHIVED_HOLD`；保持不变，不进入任何 TRAIN/DEV/TEST 集。 |
| C-D-AGL warmup 与 step1000–4000 checkpoint 序列 | `C-D-AGL-1250K-rho097-seed01`；与 step8000 完全同源、同一连续轨迹 | 日志目标为 1600 K。step0–1000 记录到约 1600–4266 K 的瞬态；后续继续至 step4000。不能把相邻 checkpoint 当作独立样本或独立来源族。 | step1000 checkpoint 记录 `QA_status=NOT_PERFORMED_BY_SAMPLER`、`STRUCTURE_PASS=false`；这些中间 checkpoint 没有各自独立的全套液态验收。末段结果只按上行的 step7000/8000 窗口评估。 | **未单独通过**；它们是一个 pilot 轨迹中的中间帧，不是独立液态候选。 | 归档 pilot；未分配 TRAIN/DEV/TEST 角色。 |
| AgL_1250K_rho0p97_seed01 | `new_liquid_Ag_1250_rho0p97_seed01` | 方案设定 1250 K；未采样 | 无结构、无 QA | 否，尚未生成，不能评估 | 仅提案；旧表拟作训练候选，须重新批准角色后才适用 |
| AgL_1250K_rho0p97_seed02 | `new_liquid_Ag_1250_rho0p97_seed02` | 方案设定 1250 K；未采样 | 无结构、无 QA | 否，尚未生成，不能评估 | 仅提案；旧表拟作训练候选，须重新批准角色后才适用 |
| AgL_1250K_rho1p00_seed01 | `new_liquid_Ag_1250_rho1p00_seed01` | 方案设定 1250 K；未采样 | 无结构、无 QA | 否，尚未生成，不能评估 | 仅提案；旧表拟作训练候选，须重新批准角色后才适用 |
| AgL_1250K_rho1p00_seed02 | `new_liquid_Ag_1250_rho1p00_seed02` | 方案设定 1250 K；未采样 | 无结构、无 QA | 否，尚未生成，不能评估 | 仅提案；旧表拟作训练候选，须重新批准角色后才适用 |
| AgL_1250K_rho1p03_seed01 | `new_liquid_Ag_1250_rho1p03_seed01` | 方案设定 1250 K；未采样 | 无结构、无 QA | 否，尚未生成，不能评估 | 仅提案；旧表拟作训练候选，须重新批准角色后才适用 |
| AgL_1250K_rho1p03_seed02 | `new_liquid_Ag_1250_rho1p03_seed02` | 方案设定 1250 K；未采样 | 无结构、无 QA | 否，尚未生成，不能评估 | 仅提案；旧表拟作训练候选，须重新批准角色后才适用 |

在本节点的 `/workspace/sampler-work` 和当前仓库检索中，没有找到第二条实际 liquid-Ag 轨迹或第二个实际候选结构。上表六个 AgL 条目来自旧规划 CSV，均明确为“proposal only; no structure generated”；C-D-AGL 的 segment archives、continuation checkpoints 和 step8000 不是彼此独立的 pilot 家族。

核对来源包括：

- `Ag64_step8000_candidate_release/MANIFEST.json`、`QA_report.json`、`thermal.csv`、`SHA256SUMS.txt`
- `/workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot/online-screen-step4000-12000/results/QA_report.json`
- `/workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot/C-D-AGL-1250K-rho097-seed01/`
- `research/mace-v12-transfer/coordination/reports/window-c/unified_potential_independent_review_20261010/liquid_ag_pilot_plan.csv`

### Ag64 step8000 的未来 DEV 候选条件

Ag64 step8000 可以保留为未来 DEV 审核候选，但当前仍为 `ARCHIVED_HOLD`，本轮不分配角色。要将来考虑 DEV，至少需要：

1. 由来源负责人解决 task ID 的 1250 K 与全程 1600 K setpoint 冲突；明确该结构代表的温度、密度与热力学状态，不能把它作为 1250 K 标签。
2. 保持它所属 C-D-AGL 来源族整体与 TRAIN 隔离；不得从同一 seed、同一初始化或其续跑中抽帧进入 TRAIN。DEV 资格由角色负责人另行审定。
3. 补齐可复核的冻结输入、sampler/势函数来源及许可、seed、环境和运行配置哈希；保留当前结构及 QA 原始哈希。候选交接包是 structure-only，不含速度、RNG 或 checkpoint；其完整来源证据需通过现有归档路径审阅。
4. 对独立平衡窗口补做统计：温度、能量、压力/应力和密度的块均值/误差；RDF 与 S(k) 尾部/峰形、配位数分布、Q6/CNA/PTM 或经验证替代方法、MSD 扩散区间；同时报告 64 原子周期盒的有限尺寸局限。
5. 最终若作为 DFT DEV 几何，使用冻结的原始坐标做统一 GPAW 单点，记录能量、原子力和应力；训练集、开发集的父族哈希审计通过后再由 C 批准角色。

## 2. 液态 Ag DFT 申请建议

### 推荐数量与来源

建议第一批申请 **6 个 64-Ag 单点结构**：1250 K 条件下三个相对密度点（经 C 冻结的参考密度的 0.97、1.00、1.03 倍），每个密度由两个独立 seed 各产生一个经 QA 筛选、时间去相关的快照。该数量沿用既有 liquid Ag 申请草案；它是申请规模，不代表结构已生成、已获准或已分配角色。

在实际结构冻结前，将六个候选限定为六个独立采样来源族；每族只选一个代表快照，避免相邻帧伪装成独立样本。若同一来源族需要多帧，它们必须整体属于同一数据角色，并按轨迹相关时间间隔去重。切勿把 C-D-AGL step8000 或其任何 checkpoint 混入这六个新 seed。

密度状态需同步一项文档证据：本任务的上游状态说明称 C 已完成密度核验，但当前 D 可见的 `unified_potential_independent_review_20261010` 仍将 1250 K 精确密度来源标为未核实。无需重下文献或参数；在任何几何体积冻结前，将 C 的已核验来源、表/方程、单位、不确定度和哈希写入可读取的项目记录，并标明 supersedes 旧表状态。

### 结构选择与来源族隔离

- 先统一并记录 sampler、版本、初始化方法、seed、密度、温度、cell、PBC、势文件 SHA-256、采样许可和轨迹索引。采样势只用于构型生成，绝不作为 DFT 标签。
- 取帧前确认温度、密度、总能量和压力/应力在多个块上的平稳性；筛选 RDF/S(k) 的液体尾部及峰形、第一 RDF 极小值处的配位数分布、Q6 和 CNA/PTM/等效局部序参量、MSD 线性扩散段；检查异常升温和晶核持续增长。
- 对周期几何进行 SHA、标准化周期指纹及 SOAP/RDF 近重复检查。一个 seed 的全部续算和派生 checkpoint 视为同一 parent/source family；同一来源族不得跨 TRAIN、DEV、TEST。
- 六个提案仅可作为未来训练候选申请；本报告不创建或改变角色。DEV 与封存 TEST 必须使用不同、独立且未进入这六个训练候选的来源族；TEST 由 C/评估负责人保管，不由 D 选帧。
- Ag64 step8000 继续 `ARCHIVED_HOLD`，不自动纳入 DEV，不作为六个提案之一，也不以同家族轨迹充当独立验证样本。

### 建议的 DFT 标签与参考设置

对每个被批准并冻结的结构，申请记录：

| 标签 | 建议记录 |
|---|---|
| Energy | GPAW 总能量及与有限电子温度/展宽一致的 force-consistent free-energy 字段；明确 eV 与能量 key 约定，可另存 eV/atom 派生值。不得与 EAM 能量混用。 |
| Forces | 所有 Ag 原子的三分量力，eV/Å；同时保存最大绝对分量和 RMS 用于结果核查。 |
| Stress | 完整 3×3 stress tensor，eV/Å³，并固定张量符号/压力约定；作为液态体积与 virial 响应的标签或至少保留为 QA 观测量。 |

沿用当前项目 DFT 基线申请：GPAW 26.7.0、ASE 3.29.0、PBE、PAW 1.2.1、PW 500 eV、电子展宽 0.10 eV、SCF 目标 1e-5。先用一个代表性固定结构做 Γ 与 2×2×2 k 点差异核验，必要时扩展网格；其余结构用审核后的统一倒格矢密度。采样势 WMH2006 Ag.eam.alloy SHA-256 为 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`，仅作历史采样 provenance，不作 DFT 标签或参考能量。

单点必须保留冻结几何哈希、参数文件哈希、PAW 身份、完整 SCF 收敛记录、能量/力/应力结果哈希和 verifier 报告。未经单独批准，不做结构优化、不变更密度/温度/势函数、不扩大到其他体系。

## 当前结论

- 实际 liquid Ag 候选：一个 C-D-AGL seed 家族；step8000 的窄范围在线筛选通过，但 1250/1600 K 元数据矛盾及完整验证缺口使其继续保持 `ARCHIVED_HOLD`。
- 实际独立 source family 数：1；可见的其他六项只是未生成方案，尚无第二条实际 liquid Ag pilot。
- DFT 申请建议：六个新 64-Ag 结构，三密度点乘两独立 seed，正式冻结 C 的密度引用和候选 QA 后再交批。
- TRAIN/DEV/TEST 角色不变；E 的 Ag–Si registry DFT 不在本报告范围内。
- 本文件是方案，不构成采样、DFT 或训练授权。
