# 是否开始受约束TSCBOP试拟合

**当前结论：`START_TRIAL_FIT=NO`。**

已有数据足以支持敏感性排序和指出TSCBOP在体相Ti/C及单个界面派生slab上明显不匹配：48原子k10×10×4的Ti/C RMSE约2.221/2.128 eV/Å；114原子TSC-only slab的Ti/Si/C RMSE约2.624/2.607/2.443 eV/Å。全列±1%探针显示A/B/λ1/λ2方向最敏感，但不同元素和体相/表面环境的偏好有冲突，未建立行参数可辨识性。

目前存在三项直接阻塞：

- Plummer–Tucker 2019原文/附录和参数约束未核实；
- bulk stress、体积/剪切响应和Si位移缺失，clean-surface数据不足；
- 没有已冻结的独立DEV/TEST family。

下一步应先收集论文来源、处理现成四个volume/shear proposal并补Si响应，再由C冻结候选数据角色和独立验证家族。若之后先做受约束小试，只能明确标记为**探索性拟合**，每次仅改变一个已核验的tied parameter block；不得写回原TSCBOP，不得调Morse抵消骨架误差，也不得声称泛化或生产可用。

本轮已停止在审计和方案阶段，等待后续授权。
