# 新窗口辅助任务：从这里开始

选择仓库 `linkengde/-`、分支 `main`，使用另一台独立云环境。先读 `coordination/README.md` 和两个任务 JSON；本窗口角色固定为 `window-b`，负责新增 DFT 数据。

## 必需软件

- Python 3.12、Git；GCC/G++、make 和对应 Python 开发头文件。
- GPAW 26.7.0（编译启用 MPI、ScaLAPACK）、gpaw-data 1.2.1、ASE 3.29.0。
- Open MPI 5.0.7、ScaLAPACK 2.2.2、OpenBLAS 0.3.29、LibXC 5.2.3。
- 本任务不需要安装 MACE 或重新下载基础模型。

先检查现有环境。仓库脚本采用当前 Debian trixie / x86-64 云环境已验证的本地依赖安装方法，先检测可用环境并跳过重装；它由本窗口实际使用的安装脚本整理而来，尚未在第二台全新机器上重复安装。

```bash
cd /workspace/-
python3 research/mace-v12-transfer/coordination/sync_tasks.py claim window-b
bash research/mace-v12-transfer/coordination/install_parallel_dft_environment.sh
```

任务领取 push 失败、任务已被别的机器领取，或者发现该任务已经完成时，暂停该任务并检查远端登记。不要继续安装后重跑同一标签。软件若缺少，可按用户既有授权补装。

## 启动分配的四个标签

```bash
cd /workspace/-
bash research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v14_parallel_acquisition/run_parallel_labels.sh
```

队列顺序：

1. `AgTi_registry_probe_v13_01`：补充 Ag–Ti 新几何检查，后续可用于 v14 训练。
2. `AgC_registry_strain_acq_v14_01`：Ag–C 新采集点。
3. `AgSi_registry_strain_acq_v14_01`：Ag–Si 新采集点。
4. `AgTi_registry_holdout_v14_01`：v14 预留盲测，禁止加入 v14 训练。

输入已生成，核对 `input_manifest.json`；不用重新运行生成器。每个点用四个 MPI 进程，线程数各为 1，保持 PW-PBE 500 eV / Gamma / 0.1 eV 方法。脚本核验输出后逐点自动提交推送；A 不必等四个全部结束才拿到前面的结果。

输出目录存在但未归档时，先检查本机进程、日志和 `.gpw`，保留现有输出并诊断恢复方法。现有启动脚本拒绝覆盖未归档目录，不能盲目重启。大 `state.gpw` 留在 `/tmp`。

此队列会逐标签自动归档并推送。若推送失败，保留当前 B 环境和本地提交；不要重算该标签，也不要修改 owner。先检查 Git 错误，在同一环境同步 `origin/main` 并解决普通 rebase/push；冲突解决后推送成功，才继续队列。A 可以并行继续自己的任务，但在标签出现于远端 `main` 前不能读取或训练使用它。

队列运行时会每 10 轮把小型进度记录推送到 `coordination/progress/window-b.json`。两边的最新进度可一起查看：

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py status --remote
```

这些新构型与现有母结构有联系，用于检验新排布和小扰动，不是完全独立的材料形貌，也不是指定温度的热平衡快照。真实周期延展界面、高温/液态环境和应变响应仍需后续专门验证。

辅助窗口只维护自己的目录和任务 JSON。当前 AgC/AgSi v13 DFT、另外两份 v14 盲测、模型训练及总日志由 A 负责。需要额外任务时先读取远端登记，等待 A 明确分配新标签。
