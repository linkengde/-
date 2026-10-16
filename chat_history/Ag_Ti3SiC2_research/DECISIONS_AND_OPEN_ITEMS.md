# 已确定决定与未解决问题（截至2026-10-16）

## 已确定

1. 物理目标：Ag局部熔池/凝固；Ti₃SiC₂固态、连通、振动及导热；不把熔坑/电弧压力作为当前主任务。
2. 已有39,486原子模型保留，但未经势验证，不能直接正式MD。
3. 旧EAM+TSCBOP+Morse存在受力风险，禁止直接用于生产模拟。
4. MACE、SevenNet、MatterSim暂只提供有限静态参考；所有候选`PRODUCTION_READY=NO`。
5. NEP/ACE开发需先严格审核DFT标签、k点收敛和独立父族；液态Ag和热态骨架/界面数据必须单独补充。
6. 本聊天生成的10个新增结构都是`UNLABELED`；1个冻结32原子基准，共11个已交接的DFT输入。
7. H窗口32原子Γ居中4组k网格的最终离线审计已完成（`cloud-h-final-kpoint-audit-20261016@2a5cd309c8be4624f98c75093c83b635aadbfa4c`）：三对相邻网格原始严格门槛**均FAIL**；事后候选0.05 eV/Å工程门槛仅后两对整体PASS，是否足够用于训练**未确定**。暂时停止这一个冻结结构的加密，1×10×8只是有限范围的候选参考，不代表k点收敛或已获训练资格。
8. 原始科学文件在科研目录，聊天内容独立存`chat_history/Ag_Ti3SiC2_research/`，不混淆文件来源。

## 开放事项与顺序

| 优先级 | 项目 | 需要的下一证据 | 当前结论 |
|---|---|---|---|
| P0 | 32原子参考k点收敛 | 最终离线审计已完成；暂停该结构加密，若恢复需独立批准新网格；候选0.05容差须另行论证 | 原始STRICT三对均FAIL；工程候选后两对PASS但适用性未定 |
| P0 | 历史1×4×2数值差异 | 追查缺失的旧k点偏移、PAW/setup、GPAW原始日志和完整参数；新旧标签严格隔离 | 新P0对旧记录Native差+0.600732 eV、力RMS差0.101692 eV/Å；未解决 |
| P1 | 纯Ti₃SiC₂骨架 | H验收V1代码/3项离线PASS但发现各MPI rank独立墙钟停止可能死锁；ChatGPT已上传V2修正为rank0广播并补5项离线测试及累计SCF/checkpoint来源日志（branch chatgpt-tsc48-segmented-runner-v2-20261016@779f9b2e），等待H复验V2。需真实4-rank写入/恢复、持久卷、墙钟、三份GPW暂存峰值+scratch预算证据。 | V3_PATCH_SUBMITTED_H_REVIEW_PENDING；REAL_RESTART_VALIDATED=NO；RESOURCE_READY=NO；DFT_EXECUTION_READY=NO；DATA_ROLE=HOLD；NOT_STARTED |
| P1 | 纯Ag相变 | 不同温度密度种子下真正液相和固液共存DFT | 不足 |
| P1 | Ag–Ti₃SiC₂界面 | 多终止面/配准/间距并定量核验法向力 | 不足 |
| P2 | 专用NEP/ACE | 严格分组训练集与盲测；小模型热力学和MD效率验证 | 尚未批准 |
| P2 | 39,486原子MD及大模型 | 全模型静态、短MD和性能/TTM校验 | 尚未获批 |

## 解释约束

- `SCF_CONVERGED`≠`KPOINT_CONVERGED`≠`POTENTIAL_VALIDATED`≠`MELTPOOL_VALIDATED`。
- 原子力RMSE需要**同一结构DFT标签**；98中心局部力数据不能替代39,486原子全体结构。
- 不得将Ag32固态随机位移当作液态Ag，不得将同父界面不同gap看成独立盲测。
- 任何云端计算启动/结束，只能以真实执行日志、结果目录和校验记录声称。
