# MACE V17 单主熔池续算分支（2026-10-22）

从原MACE V17紧凑试算的已提交检查点恢复；原始模型、权重、轨迹、检查点均保存在上级 `mace_v17_single_pool_bicontinuous_20261021/`，此目录是隔离副本。

- 模型：766原子紧凑双连续试算模型（不是37,711原子目标模型）。
- 势：MACE V17，SHA-256由 `run_compact_mace.py` 启动时强制校验。
- 续算起点：累计加热0.81 ps；`HEATING_RESUME_CHECKPOINT.extxyz` 与 `HEATING_RESUME_STATE.json` 成对保存，含位置、速度和Langevin随机数状态。
- 设置：0.5 fs，Ag局部1900 K Langevin源，Ti₃SiC₂热沉；每100步保存一次检查点，每10步统计。
- 热源区外液态候选团簇继续记录，但不用于停机，也不作为主熔池证据；熔化判据只检查受热Ag区的持续连通、q6无序、位移和候选液态团簇。
- 保留物理加热上限4 ps；没有墙钟时间上限。严重力/短接触仍触发安全停止。若多指标主熔化判据持续满足，将自动转入冷却。
- 输出只写本目录，不覆盖上级原始结果。

OVITO轨迹：`HEAT_COOL_OVITO_RESUMED.extxyz`。运行：在此目录使用 `/tmp/mace-v17-cpu-env/bin/python -u continue_heating_chunked.py`。
