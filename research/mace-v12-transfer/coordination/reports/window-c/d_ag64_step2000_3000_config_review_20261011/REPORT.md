# Ag64 step2000→3000 候选签署审核

结论：**REJECTED**。签署配置未生成，`SIGNED_CONFIG_SHA256=NOT_ISSUED`；`execution_authorized=false`。没有启动真实 MD、DFT 或训练。

审核来源：分支 ag64-step2000-extension-dev，提交 538f23e6fc2b3b3461693085487e065325229e5e。实际配置 SHA 匹配用户给定 1aa03f05ed2a3662985269c2b0148542d79202f771aff690d67894e25f7c1f96。

## 六项审核

1. step2000 归档、state 和 state_arrays 实际哈希均匹配。C 在 ASE 3.29.0 下读取完整归档，恢复64 Ag、坐标/动量与 PCG64 状态，PASS；未附加计算器或调用积分。
2. 势 SHA pin 与冻结 Ag.eam.alloy 一致。D 本机势字节仍须运行前实时核验。
3. adapter/core/backend 实际 SHA 均匹配候选记录。
4. ENVIRONMENT_SHA 字段已记录，但没有可复算的环境 manifest。候选引用的 fff2d118 提交实际为 ASE 缺失报告，明确写 dry-run NOT_RUN，与候选 PASS 声明冲突。C 的独立恢复检查仅补充 C 本机证据，不能替代 D 实际环境身份。
5. start=2000、end=3000、max_added=1000 字段正确；ASE step 方法在3000之前拒绝第3001步，独立边界检查 PASS。已有6项 mock测试复跑6/6通过。它们不覆盖真实积分后的中断恢复或 backend v2 全链路。
6. 暂不能签署：需要补齐环境身份、backend v2版本对应的非科学性远端往返证据以及超时恢复能力。

## 用户追加的墙钟预算调整

按用户“至少30分钟”的要求，新未签署候选采用3600秒（1小时）资源上限。只在独立 resource_budget_revision 目录保存新代码副本与候选：adapter的预算校验、core预算默认值/绑定和JSON两处预算同步更新；原候选、原执行器、checkpoint、势和全部物理参数保持原样。所有新增代码 SHA 已重新计算并写入新版候选。变化仅是预算，执行范围始终2000→3000、最多1000步。

3600秒是保守资源预算，不能保证完成。历史约605秒/300步包含归档开销，不能直接当作纯积分性能；新归档频率改变后需实际性能证据。

## 超时恢复阻塞

新core使用1000步归档周期，而source sampling字典仍记录 checkpoint_every=100。advance_segment_with_outputs 在尚未到第一个1000步归档点时超时，直接抛 SEGMENT_WALLTIME_BUDGET_EXCEEDED，没有终止检查点保存逻辑，也无SIGINT/SIGTERM恢复处理。独立无原子、无积分的计数器故障注入证明超时出口未产生文件，见 TIMEOUT_STATE_MACHINE_PROBE.json。旧step2000不会被覆盖，但临时工作区中新积分状态可能丢失。加长预算不能消除这一故障路径。

此外，新 backend v2（SHA 0a8d8b71ee938e768f47bc9e819c5953f7f7bb0aa6415f248916b4a300a3928f）未提供版本对应的模拟输出上传、全新fetch、恢复及状态比较报告。旧backend通过记录不能自动证明新SHA版本通过。

## 重审所需证据

- 可复算的ENVIRONMENT_SHA原始payload，D实际ASE/NumPy/Python/runtime信息与真实只读恢复报告，修正候选引用错误。
- 与最终代码SHA一致的backend v2合成upload/verify/restore报告；覆盖远端失败、成员SHA和positions/velocities/RNG恢复。
- 超时/中断后原子化保留并可靠备份最后完整步的方案与故障测试，明确100与1000检查点间隔角色；不得静默修改物理采样节奏。
- 3600秒修订资源版本重新冻结SHA；完整审核通过前不生成执行授权配置，不运行MD。
