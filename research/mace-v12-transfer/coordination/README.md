# A/B 窗口协作：当前分工

GitHub `main` 同步输入、任务登记、紧凑结果和日志；实时进程、Python 环境与 `/tmp` 检查点仍在各自云机器。A 在 B 当前任务完成前把下一项不冲突的任务写入 `tasks/window-b.json`。B 会在当前任务发布后同步 `main` 并接续下一项，不需要用户逐条转发命令。仓库提交不能唤醒已经结束的 B 会话；结束后只需恢复该环境。

## 当前状态与任务

| 窗口 | 当前工作 | 不负责的工作 |
|---|---|---|
| A (`40f189e9-44af-4d33-8992-e3b7b44dbf3b`) | V17 development-validation：AgC 已收敛且归档核验通过；AgSi 四核运行中，AgTi 排队。随后审阅 B 的训练来源审计和三项开发验证结果，再构建、训练、评估 V17。 | 六份 V17 withheld-test 标签继续封存，直到所选 checkpoint 和选择记录冻结。 |
| B (`c5035b48-f43b-4da0-b8e4-e2862817f86a`) | 正在审计 35 个 V16 训练帧和 8 个 V17 配对训练诊断标签的来源/能量约定；随后自动接续 V17 分割泄漏几何审计。 | 不读取 A 的开发验证输出或六份 withheld-test 标签；不擅自运行 DFT、训练/推理、MD/TTM。 |

A 的 DFT 和 B 的 CPU 分析在独立云环境中并行。当前 v14 三份盲测的力向量 RMSE 都高于 0.05 eV/Å，因此 v14 不用于长时间 MD/TTM。精确进度见 `coordination/tasks/window-*.json` 和 `coordination/progress/window-*.json`。

## B 当前任务

无需为队列中的每个任务单独给 B 发命令。A 维护任务卡和后续队列；活跃的 B 会话按 `B_CONTINUOUS_WORK.md` 自动继续。若 B 会话已结束，恢复该云环境后它再同步仓库。若 owner 校验失败或身份不符，B 停止并回报；不要手动改 owner。

## 同步与归档

- 开始任务前先 `git fetch origin main` 并快进同步。B 已登记任务只能由同一云实例继续；身份不符时不接管、不重复计算。
- 每次 DFT 标签须先通过输入哈希、结构/ID、能量/受力有限性、摘要和 GPAW 日志核验，再发布。进度 watcher 每 10 次 SCF 更新一次，MACE watcher 每 4 个 epoch 更新一次。
- A 拥有主工作日志和全局 `file_manifest.json`；B 只发布自己的 `coordination/reports/window-b/`、B task/progress 和已授权的计算产物。
- 不上传大型 `state.gpw`。推送/rebase 失败时保留原 owner 环境、检查点和结果；先正常同步 `main`，不得 force push 或重算已完成标签。
- v15/v16 必须按前一版误差补数据，并各自保留全新的盲测。现有派生小簇构型不能证明扩展周期界面或高温/液态 Ag 的迁移适用性。模型通过逐界面验证和结构稳定性关卡前，不做生产长 MD/TTM。

远端进度查看：

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py status --remote
```

## Standing B work queue

See [B_CONTINUOUS_WORK.md](B_CONTINUOUS_WORK.md). A replenishes the queue in
tasks/window-b.json before B's current task ends when practical; only one B task
runs at a time. B publishes stages, syncs main and continues assigned work while
its session remains active. Git alone cannot restart an ended agent session.
