# V19 受控训练启动门

状态：**未获训练授权；本次未启动训练。** 推荐起点仍是原始 MACE-MP-0b3-medium；冻结 V18 只作固定对照，不从 V18 权重继续微调。该门槛沿用 V19 修订方案，并增加液相能量口径要求。

## 启动前必须满足

1. **训练角色与来源合法。** 每个 frame 有可审查 owner、源结构/DFT input/结果 SHA、确切 SCF+verifier 状态、PBE/PAW/PW/smearing/kmesh、forces 与 energy key。当前 15 Ag/10 Ti₃SiC₂ 仍是数值控制，不进入建集；只有另行授权明确批准逐条角色变更后才可用。V18 原 30 train 可以作为受控 R30 输入复用；V18 3 dev 继续标为历史已见 development。
2. **能量—力一致。** 同一拟合包 energy 和 forces 同源、eV/eV·Å⁻¹、固定原子顺序/ID/cell/PBC。V18 recipe reproduction 全部用同一个原生 `REF_energy`/`REF_forces` key。生产 thermal fit 使用 force-consistent `free_energy` 与该次有限展宽力；若 V18 训练帧无法从已验证归档重建同口径 free-energy，则禁止直接将其与 thermal frames 混训练，需单独授权重整/排除策略。禁止 native/free energy 混训。
3. **parent family 拆分先冻结。** 按 source geometry/CIF hash、cell/PBC、原子映射、trajectory seed、interface registry/termination 和全部派生 history 建 manifest；同一父结构的位移、k-mesh、热帧、gap scan、slab 数值控制整组只留在一个 split。已评分 22 个 C 静态结构、六个 residual 接触族、V18 三个历史 dev 父族永久不可作为独立测试。
4. **开发集独立可见，测试集封存。** 训练启动前必须有新 parent-disjoint dev 数据可供 early stopping/hyperparameter 选择；测试结构/标签由独立 evaluator 持有，manifest hash 预先冻结，C 不读。若只能用同父族开发数据或封存数据来选择模型，停止，不训练。
5. **骨架无灾难性退化。** 提案门槛：新臂在同一 parent-disjoint Ti₃SiC₂ dev 上向量力 RMSE ≤0.05 eV/Å，且不得超过同构型 base MACE 的 1.25 倍；Ti、C 的 z RMSE 和 max per-atom error 单列，禁止由 Si/Ag 或正负层误差抵消。对称位移恢复斜率相对 DFT 误差暂拟 ≤20%，正式值须在独立 dev/数据覆盖核实后冻结。若任一骨架门失守，即停止/否决训练臂，不用增加 epoch 掩盖。
6. **Ag-X 接触改善且无跨域交换。** 对 Ag-Ti、Ag-Si、Ag-C 分别报告向量 RMSE、元素 xyz、Ag-X 分离投影、Ag 层 z 广义力和最大逐原子误差。提案门槛：独立 dev 每个接触向量 RMSE ≤0.05 eV/Å、分离力绝对误差 ≤0.10 eV/Å；同时确认相对于 base MACE 有实质改善，且骨架门仍通过。旧混合 EAM+Tersoff+Morse 的误差单独记录，不并进 MACE 误差。
7. **训练实验可归因。** 先用基础模型独立初始化，保留 R30、冻结低层 F30、骨架回放 R30+S 的单因素对照；匹配 seed、epochs/update budget、训练器版本、参数/梯度 mask、loss 和输出目录。当前 R30+S 在数据角色门获批前不能组装；F30 在 MACE freeze API/parameter-hash 测试通过前不能跑。所有臂独立目录，不覆盖 V18。
8. **相变数据单独接入。** 六帧 Ag 液体只作 1250 K 近熔候选，需液态判据/多 seed/密度来源与 DFT mesh 门通过。未来至少还需其它温度/密度液体、过冷液体、多个固液共存取向/尺寸、真实热界面、高温骨架和近程排斥标签。报告模型训练 force score 通过不能替代熔点、密度、潜热、扩散、RDF、界面迁移和冷却速率检验。

## 运行期停止规则

- 任一训练标签 SCF/hash/energy key/atom order/parent family/PBC 不清：移出候选包并停止 preflight。
- dev 骨架力超过上述门或明显劣于基础模型：立即判该臂不合格，检查数据/损失归因，不以训练时长补救。
- interface 改善以 Ti/C 骨架恶化换取，或某个 Ag-X 类别仍超过门：整体不通过，不能只报告 pooled RMSE。
- 测试标签被查看、dev/test parent 交叉、训练结果将覆盖 V18/A/B 文件：停止流程并保留证据。
- MACE API 不支持设计中的 freezing 或 update budget 不能匹配：冻结臂标为不可执行，先做独立实现审计，不塞入未经核验的 CLI 参数。
- 未完成液态/共存物性验证时，模型最高称“固相/界面开发候选”，不能称为熔池—凝固或熔坑生产势。

上述阈值是待批准的研究门槛，不是已通过的结果，也不自动授权训练。
