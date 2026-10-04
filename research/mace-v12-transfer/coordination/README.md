# 两个云窗口的分工与同步

GitHub `main` 是交接依据。它同步输入、任务登记、结果和日志；实时进程、虚拟环境、`/tmp` 检查点仍各自在对应云机器上。

## 固定分工

| 任务 | 负责人 | 输入与结果目录 | 工作 |
|---|---|---|---|
| `window-a` | 当前窗口 | `periodic_interface_v4/pbe_interface_v13_holdouts/`；`pbe_interface_v14_main_holdouts/`；模型目录 | 完成已启动 AgC/AgSi DFT，再计算 v14 AgC/AgSi 预留盲测；统一构建、训练、评估 v14，必要时继续 v15/v16；维护主日志与总清单 |
| `window-b` | 新的独立云环境 | `periodic_interface_v4/pbe_interface_v14_parallel_acquisition/` | 计算三份 AgTi/AgC/AgSi 新采集标签，再计算一份 AgTi 的 v14 预留盲测；逐点核验并提交小型结果 |

窗口 B 不重跑窗口 A 的两份 DFT，不重训 v12/v13。窗口 A 不启动窗口 B 已领取的标签。窗口 B 只更新自己的结果目录和 `tasks/window-b.json`；主工作日志、总清单和模型训练由 A 汇总，避免修改相同文件。

## 开工与任务领取

1. `git fetch origin main`，在保留本地改动的前提下同步最新 main。
2. 阅读本文件、`tasks/window-a.json`、`tasks/window-b.json`、主工作日志及最新模型评估。
3. 新窗口运行 `python3 research/mace-v12-transfer/coordination/sync_tasks.py claim window-b`。普通 Git push 成功后才启动计算；并发领取时，第二次 push 会被拒绝，拒绝的窗口不能计算。
4. 已有 owner 的任务仅由同一云机器继续。过期时间戳不证明作业已经停止；不得自动接管或重复启动。

本窗口任务已登记。`owner_instance` 是保存在 `/workspace` 的随机机器标识，不是凭据。共用同一机器的窗口会看到相同标识，不能领取需要另一台机器的 B 任务。

## 结果回传与使用

窗口 B 的启动脚本在每个标签收敛后核验输入、结构、有限标签和日志；写入 `verification.json`，随后自动提交并推送自己的结果和任务进度。`.gpw` 不入库。推送失败或 rebase 冲突时保留结果并暂停队列，确认同步成功后再继续。

读取另一窗口进度用：

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py status --remote
```

窗口 A 在构建下一版数据前拉取 main 并重新核验 B 的结果；此前的盲测点转入下一版训练后，不再作为该版独立盲测。每版仍有未训练的验证结构。v15/v16 的补点和训练策略根据上一版误差决定，不预先重复同一套训练。

两边各承担四份 DFT 单点（A 包括当前已运行的两份），让计算量较均衡。A 的模型工作在 DFT 作业之间安排；同一台四核机器上不同时运行四进程 DFT 与 MACE 训练。

若需要重新分配任务，由 A 在确认原 owner 的作业停止或已完成后更新登记并推送。不得用 force push 抢占任务。
