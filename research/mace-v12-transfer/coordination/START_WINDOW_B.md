# 窗口 B 的协作任务

B 是独立云窗口。GitHub 只同步任务和结果；A 每次安排 B 开始新阶段时，仍需把对应命令发到 B 窗口。四份 DFT 标签已全部完成，不要重算。

## 现在执行：定位 v13 受力残差

在 B 窗口执行下面命令。脚本会核对/准备 MACE CPU 环境、检查五份已验证 acquisition 标签，生成每个结构的逐原子误差、界面附近/远端骨架误差及最高误差原子的邻居，并将报告推回仓库。报告用于决定后续补点，不改变 v14 数据和盲测。

```bash
cd /workspace/-
git fetch origin main
git merge --ff-only origin/main
bash research/mace-v12-transfer/coordination/reports/window-b/run_v13_force_localization.sh
```

## 稍后执行：v14 训练和评估

等 A 在 `main` 发布 `mace_periodic_v14_interface_energy/data/dataset_manifest.json` 后，再把下面命令发给 B 窗口。不要在数据清单出现前训练。入口会重新核对 31/2/3 划分、4/4 能量组成秩、数据哈希和三份 v14 盲测核验/隔离；随后只训练一次并运行冻结测试与三份盲测评估，实时同步 epoch。

```bash
cd /workspace/-
git fetch origin main
git merge --ff-only origin/main
bash research/mace-v12-transfer/coordination/reports/window-b/run_v14_training_and_evaluation.sh
```

B 任务仍归 `c5035b48-f43b-4da0-b8e4-e2862817f86a` 所有；身份校验失败时停止，不改 owner。B 不重训 v12/v13、不修改主工作日志或总 `file_manifest.json`，也不运行生产长 MD/TTM。

## 已完成的 DFT 任务（历史记录）

B 的四个标签 AgTi probe、AgC strain、AgSi strain 和 AgTi v14 holdout 已全部收敛、核验并推送。历史计算输入/设置仍保存在 `periodic_interface_v4/pbe_interface_v14_parallel_acquisition/`；不要再运行旧的 `install_parallel_dft_environment.sh` 或 `run_parallel_labels.sh`。

这些新构型与现有母结构有联系，用于检验新排布和小扰动，不是完全独立的材料形貌，也不是指定温度的热平衡快照。真实周期延展界面、高温/液态环境和应变响应仍需后续专门验证。

B 现负责 v13 受力误差定位和 v14 训练/评估；A 负责剩余 DFT、v14 数据集构建、主日志和总清单。开始每个 B 阶段都要由 A 给出一条具体命令。
