# 势函数开发下一步决策

## 决策

- `CURRENT_TSCBOP_FOR_PRODUCTION=NO`。实现和原子映射通过，但 TSCBOP-only 在 135 原子组件力比较中 RMSE 为 **2.5624 eV/Å**，在独立 48 原子体相也有约 **1.99 eV/Å**。不能把当前势作为已验证 Ti₃SiC₂ 骨架势。
- `ERROR_SOURCE= BULK_BACKBONE + SURFACE/INTERFACE/RELAXATION`。bulk 配位样核心仍有残差；终止 Si、界面相邻 Ti 及局部面内重排区更差；底部欠配位 Ti 反而较低。135 原子模型不能完整表征厚体相，误差既不能归为纯 bulk，也不能归为纯 surface。
- `PRIMARY_SCIENTIFIC_ROUTE=MACE`。目标是同一 Ag–Ti₃SiC₂ 多元素势；在足量 PBE/PAW bulk、strain、surface、registry/gap 与温度样本下，DFT 训练的多体势最灵活。MACE 需要独立 TRAIN/DEV/TEST 设计及 LAMMPS 接口/力能量一致性验收。
- `TSCBOP_REFIT=HOLD`。先取得原始论文与参数来源、完成已有四个体积/剪切提案的 DFT 单点，再决定参数化形式是否值得保留。禁止通过 Morse cross 参数抵消 Ti/C 骨架误差。
- `ALTERNATIVE_CLASSICAL_POTENTIAL=NO_VERIFIED_CANDIDATE`。本地档案没有可凭来源与测试推荐的第二套 Ti–Si–C 参数；不以势函数类型名替代验证。
- `875K_ATOM_PRODUCTION_ON_CURRENT_ENV=NO`。当前快照是 4 CPU quota、32 GiB、无 GPU。短程经典势或许能装下，但吞吐未实测；当前 CPU-only MACE 也没有 875k 的并行、内存和速度基准。应在有 GPU/多核持久化节点上分级 benchmark，再考虑生产任务。

## 按顺序执行的工作

1. 保持原始 TSCBOP/EAM/Morse 文件和已完成 DFT 数据不变；先解决 TSCBOP 论文/参数来源不可核验问题。
2. 在 C 冻结角色与正式参数后，优先完成已存在的四个 bulk 体积/剪切 DFT proposal。若 E/F 资源不足则转到更合适的持久 DFT 节点，不降低精度。
3. 基于已有 bulk baseline、4 个 ±0.02 Å 位移和新增体积/剪切响应，对 TSCBOP 做一次可行性决策；不是自动拟合。若骨架残差仍明显，停止 A 路线并转 MACE/新参数集。
4. 用已归档 114 原子 slab 做初步 surface 诊断。只有在 clean surface 物理是目标且有已核验结构来源时，再补独立 slab DFT。
5. 骨架方案通过后，由 C 冻结多注册、多个 gap、近平衡扰动与分离参考，开展 Ag cross screening。保留当前 2.5624、1.3513、1.5178 eV/Å 的既有组件结果为基线，不与后续结果混用或覆盖。
6. 在约 875,000 原子模型上，先在 10k 和 100k 原子规模 benchmark CPU/GPU、MPI、RAM、吞吐、checkpoint 和 LAMMPS pair-style/TTM 接口，再决定生产硬件。

## 最终问答

1. **TSCBOP 是否适合当前模型？** 以当前哈希和已测试精度，不适合直接做预测性生产模拟；实现正确不改变这个判断。
2. **误差主要来自哪里？** 体相 Ti/C 骨架偏差已经存在，界面相邻 Ti、Si 终止面和局部重排区会进一步放大；单一 135 原子 slab 不足以作唯一归因。
3. **应优先重标定、替代还是机器学习？** 科学目标上优先建立 DFT 训练的 MACE 路线；先用四个 bulk response 点做数据缺口补齐。A 仅作经论文核验后的低成本对照，B 目前没有有证据的候选。
4. **Ag cross 最少需要什么？** 除保留当前基准，先增加 2 个非等价注册 × 3 个 gap（6 个）、2 个确定性近平衡位移、1 个大间距分离参考；C 冻结数据角色后执行。若要持续做分项归因，在代表点加同几何 phase-only 参考。
5. **哪条路线符合现有资源？** 当前 4 核/32 GiB/CPU-only 适合复用小规模 DFT 结果与准备/测试，不足以认可 875k 生产吞吐。经典势成本最低但该 TSCBOP 精度失败；MACE 是科学首选，生产需更强 GPU/多核资源并先做规模基准。

**本轮状态：** 只读取现有归档并生成诊断/计划文件；未改原始数据、势文件、结构或角色，未运行 DFT、LAMMPS、MD、TTM、MACE 训练/推理，也未提交或推送 Git。
