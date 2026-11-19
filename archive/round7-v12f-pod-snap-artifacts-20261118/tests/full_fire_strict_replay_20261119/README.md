# 347,464 原子 FIRE 严格 lost-atoms 复核（2026-11-19）

使用原始 `start_wrapped.data`、POD v12f、`boundary p p m`、5 个 MPI rank、z=0–20 Å 固定组、FIRE `dmax=0.02 Å`、`ftol=0.5 eV/Å`，并保留 `thermo_modify lost error`。迭代上限为 62；每步输出原子 ID 和坐标。

## 结果

- 严格 `lost error` 重放在第 62 步正常结束，退出码 0；没有 `Lost atoms` 错误。
- 第 60、61、62 帧各有 347,464 行、347,464 个唯一 ID；没有重复、缺失或非有限坐标。z 方向没有越界。x/y 有少数坐标比主盒边界略超 0.032 Å 以内；x/y 是周期方向。
- 第 60 步与此前 `lost ignore` 重放的逐原子坐标完全相同（最大绝对坐标差 0 Å），且势能 −2,237,746.355850081 eV、Fmax 49.701958859 eV/Å 与原报错运行打印值一致。
- FIRE 因达到 62 次迭代上限而停止，未收敛。

因此原始报告没有复现成原子真实逃逸；确切的 lost-atoms 内部触发原因仍未确定。间歇性 MPI/周期边界原子迁移或计数问题是可能解释，但目前没有证据确认具体机制。不要把这个未收敛 FIRE 终态用于 MD。

`strict_replay_frames_60_62.dump.gz` 保存了 60–62 步全部原子 ID、类型和坐标。完整 0–62 逐步 dump（675 MiB）留在原工作区；SHA-256 记录于 `strict_replay_audit.json`。
