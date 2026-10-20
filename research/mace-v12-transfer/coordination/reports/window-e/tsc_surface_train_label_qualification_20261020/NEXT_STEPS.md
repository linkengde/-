# 后续数据步骤

1. 保持四个 COD 9009647 表面点为 `TRAIN_FORCE_LABEL_CANDIDATE`，将已核验的质量信息和数值 mixer 差异登记到独立用途映射；不拆分参考点/位移点，不修改原始角色文件。
2. DEV：先查找是否有真正独立且已有 DFT 标签的 Ti₃SiC₂ 父族。E 当前任务清单中 DEV 参考点及其响应点均为空并处于 HOLD。若没有可复用来源，再选取有可追溯来源的独立父结构，冻结结构/输入 SHA 后按已定两点预算完成 DEV 参考与响应单点。
3. TEST：按既有任务清单由独立 evaluator 负责两点。E 不接触 TEST 结构身份或标签；确认 evaluator 已有独立父族结果即可，不将 COD 9009647 的任何点转作 TEST。
4. DEV/TEST 父族和标签来源可核验后，冻结统一的 TRAIN/DEV/TEST manifest。训练仅使用获准 TRAIN 标签；DEV 用于选择，封存 TEST 只做最终评价。当前不启动拟合/训练，不修改 TSCBOP、Ag EAM 或 Morse。

本轮只完成标签完整性与父族审计。没有生成新结构，也没有启动 DFT、MD 或势函数拟合。
