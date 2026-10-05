# 两个云窗口的分工与同步

GitHub `main` 是交接依据。它同步输入、任务登记、结果和日志；实时进程、虚拟环境、`/tmp` 检查点仍各自在对应云机器上。

## 固定分工

| 任务 | 负责人 | 输入与结果目录 | 工作 |
|---|---|---|---|
| `window-a` | 当前窗口 | `periodic_interface_v4/pbe_interface_v13_holdouts/`；`pbe_interface_v14_main_holdouts/`；模型目录 | 完成 A 的 v14 AgC/AgSi 盲测并构建、核验 v14 数据集；审查 B 的 v14 训练/评估，维护主日志与总清单 |
| `window-b` | 独立云环境 | `coordination/reports/window-b/`；`mace_periodic_v14_interface_energy/` | DFT 标签已完成；先定位 v13 逐原子受力残差，再用独立环境训练并评估 v14 |

窗口 B 不重跑任何已完成的 DFT，也不重训 v12/v13。窗口 A 不启动窗口 B 已领取的标签。当前 A 继续完成 DFT、构建/核验 v14 数据集并维护主工作日志和总清单；B 独立定位 v13 受力残差，并在 A 发布完备数据集后负责唯一一次 v14 训练和盲测评估。B 只更新自己的分析报告、v14 训练产物和 `tasks/window-b.json`；A 审查并维护主日志及全局清单，避免两边修改同一文件。

## 开工与任务领取

1. `git fetch origin main`，在保留本地改动的前提下同步最新 main。
2. 阅读本文件、`tasks/window-a.json`、`tasks/window-b.json`、主工作日志及最新模型评估。
3. 新窗口运行 `python3 research/mace-v12-transfer/coordination/sync_tasks.py claim window-b`。普通 Git push 成功后才启动计算；并发领取时，第二次 push 会被拒绝，拒绝的窗口不能计算。
4. 已有 owner 的任务仅由同一云机器继续。过期时间戳不证明作业已经停止；不得自动接管或重复启动。

本窗口任务已登记。`owner_instance` 是保存在 `/workspace` 的随机机器标识，不是凭据。共用同一机器的窗口会看到相同标识，不能领取需要另一台机器的 B 任务。

## 结果回传与使用

窗口 B 的启动脚本在每个标签收敛后核验输入、结构、有限标签和日志；写入 `verification.json`，随后自动提交并推送自己的结果和任务进度。`.gpw` 不入库。推送失败或 rebase 冲突时保留结果并暂停队列，确认同步成功后再继续。

两边不需要等对方整批结束：谁先完成，就先把自己的已核验标签推送到 `main`。A 的 v13 队列在两份标签都归档后自动运行 `pbe_interface_v13_holdouts/publish_independent_holdouts.py`；该入口核对两个 `verification.json`、归档清单及输入/摘要/输出哈希，更新工作日志和总哈希清单，并只提交 A 的 v13 结果目录及相应状态文件。它拒绝 `.gpw` 和 100 MB 及以上文件。B 的队列仍按标签逐个提交。

## 新的 B 窗口协作任务

B 的四份 DFT 已全部完成。B 现在接手 v13 逐原子受力误差定位和 MACE 环境预检；A 发布通过门槛的 v14 数据集后，B 再用独立环境完成唯一一次 v14 训练和盲测评估。每次切换阶段都必须给 B 窗口一条明确命令；GitHub 任务/进度记录只同步状态，不会自动启动另一窗口。命令和硬门槛见 `coordination/START_WINDOW_B.md`。

## 新的 B 窗口协作任务

B 的四份 DFT 已全部完成。B 现在接手 v13 逐原子受力误差定位和 MACE 环境预检；A 发布通过门槛的 v14 数据集后，B 再用独立环境完成唯一一次 v14 训练和盲测评估。每次切换阶段都必须给 B 窗口一条明确命令；GitHub 任务/进度记录只同步状态，不会自动启动另一窗口。命令和硬门槛见 `coordination/START_WINDOW_B.md`。

运行中的小进度由队列旁的 watcher 每增加 10 个 SCF 迭代推送一次；MACE 训练按 4 个 epoch 更新一次。记录任务名、状态、轮数、耗时和时间戳。两个窗口分别写 `coordination/progress/window-a.json` 与 `window-b.json`，所以并发时不会改同一份进度文件；`status --remote` 一次显示两边任务和最近进度。每个文件保留最近 100 条事件。实时日志和任务状态是可变文件，不进入总哈希清单；完整 GPAW 日志和 `.gpw` 检查点仍留在计算机本地，归档后才推送紧凑结果。

另一窗口只能使用已经出现在 `main` 且核验为 PASS 的结果。若计算完成但推送失败，计算 owner 保留原云环境、输出和本地提交；不重算、不 force push，也不让另一窗口接管。先按失败信息同步 `origin/main` 并解决普通 rebase/push；有 rebase 冲突时在原 owner 环境解决后继续推送。未推送标签在远端仍保持未完成，依赖它的评估/数据集先等待；另一窗口可以继续不依赖该标签的已分配工作。若原环境确实不可恢复，由 A 先确认进程停止、检查点/输出是否可恢复，再更新 owner 并推送后重新分配。

读取另一窗口进度用：

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py status --remote
```

窗口 A 在构建下一版数据前拉取 main 并重新核验 B 的结果；此前的盲测点转入下一版训练后，不再作为该版独立盲测。每版仍有未训练的验证结构。v15/v16 的补点和训练策略根据上一版误差决定，不预先重复同一套训练。

两边已各承担四份 DFT 单点（A 包括当前 AgSi 盲测）。接下来转为跨类型并行：A 完成该盲测和数据集审计，B 做 v13 误差定位及 v14 训练/评估。不得在 A 的四进程 DFT 环境同时开 MACE 训练；B 使用自己的独立云环境。

若需要重新分配任务，由 A 在确认原 owner 的作业停止或已完成后更新登记并推送。不得用 force push 抢占任务。
