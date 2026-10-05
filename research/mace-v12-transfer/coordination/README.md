# A/B 窗口协作：当前分工

GitHub `main` 同步输入、任务登记、紧凑结果和日志；实时进程、Python 环境与 `/tmp` 检查点仍在各自云机器。切换 B 的阶段时，A 必须把 `START_WINDOW_B.md` 中的具体命令发给 B；仓库更新不会自动启动另一个窗口。

## 当前状态与任务

| 窗口 | 当前工作 | 不负责的工作 |
|---|---|---|
| A (`40f189e9-44af-4d33-8992-e3b7b44dbf3b`) | v14 `SCREEN_FAIL`。V15 Ag-C、Ag-Si 残差标签已核验；Ag-Ti 残差 DFT 正在四核运行，最新日志到 SCF 第 20 轮。完成后由 A 依次计算三个冻结的 V15 registry 盲测点，再独立构建、训练和评估 V15。 | 不计算 B 的审计任务；三个 V15 盲测不进入训练/验证/模型选择。 |
| B (`c5035b48-f43b-4da0-b8e4-e2862817f86a`) | 几何审查已完成。当前分配只读 V15 数据完整性审计：标签来源、切分泄漏、旧 LCAO 力标签一致性及数据构建验收条件。 | 不运行 DFT/MACE/MD/TTM；不修改训练数据、脚本、标签或盲测输入。 |

A 的 DFT 和 B 的 CPU 分析在独立云环境中并行。当前 v14 三份盲测的力向量 RMSE 都高于 0.05 eV/Å，因此 v14 不用于长时间 MD/TTM。精确进度见 `coordination/tasks/window-*.json` 和 `coordination/progress/window-*.json`。

## B 当前命令

把 `coordination/START_WINDOW_B.md` 中的同步、身份检查和进度命令发给 B。B 只发布只读审计报告和哈希清单。若 owner 校验失败或身份不符，B 停止并回报；不要手动改 owner。

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
