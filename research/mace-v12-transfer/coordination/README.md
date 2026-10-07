# A/B 窗口协作：当前分工

GitHub `main` 同步输入、任务登记、紧凑结果和日志；实时进程、Python 环境与 `/tmp` 检查点仍在各自云机器。A 在 B 当前任务完成前把下一项不冲突的任务写入 `tasks/window-b.json`。B 会在当前任务发布后同步 `main` 并接续下一项，不需要用户逐条转发命令。仓库提交不能唤醒已经结束的 B 会话；结束后只需恢复该环境。

## 当前状态与任务

| 窗口 | 当前工作 | 不负责的工作 |
|---|---|---|
| A (`40f189e9-44af-4d33-8992-e3b7b44dbf3b`) | V17 provisional 43/3 data rebuilt with complete B archive manifest; preflight PASS. Run02 completed80 but exported epoch72 and is ineligible under the frozen epoch79 rule; corrected run03 starts with all checkpoints retained and an epoch79 guard. | Six V17 withheld-test labels stay sealed. |
| B (`c5035b48-f43b-4da0-b8e4-e2862817f86a`) | Archive reconciliation, approved frame0 DFT (50 SCF; verification PASS) and independent 43-row audit are complete. Next: isolated certified30 builder draft; then V9–V17 review after A’s valid fixed-epoch score. | No development reference values or withheld outputs; no edits to A’s dataset or unauthorized calculations. |

A 的 DFT 和 B 的 CPU 分析在独立云环境中并行。V16 冻结模型在三个 V17 开发结构上的力向量 RMSE 为 AgC 0.1376、AgSi 0.0729、AgTi 0.1015 eV/Å；分离力绝对误差为 0.2044、0.2035、0.1856 eV/Å，均超过暂定门槛。开发评分不是独立最终盲测，也不支持生产 MD/TTM。精确进度见 `coordination/tasks/window-*.json` 和 `coordination/progress/window-*.json`。

新开发结构的逐原子误差定位显示，AgC 与 AgTi 的 4.5 Å 接触壳误差高于壳外；AgSi 的误差较分散，但标记 Si 与邻近 Ag 存在较大的横向分量。A 已将逐原子证据发布在 V17 acquisition 目录。该结果支持先验证已完成的 V17 配对训练数据覆盖效果；不能据此宣称已找到唯一因果或直接调整损失权重。

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
