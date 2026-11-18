# Round7 POD/SNAP v12f–v12j 归档

最新状态见 [ROUND7_CURRENT_STATUS_20261118.md](ROUND7_CURRENT_STATUS_20261118.md) 和 [v12j 复核记录](v12j/ROUND7_STATUS_V12J_20261118.md)。v12j 保持 128/12 训练验证划分与 Si01 留出，加入 1 个收敛的 36 原子 Ti–C DFT 局部标签。12 帧验证 RMSE 为 POD E/F 0.077857 eV/atom、0.590685 eV/Å；SNAP E/F 0.063033 eV/atom、0.722126 eV/Å。

完整 458,248 原子 POD 冷态 10 fs 温度从 300 K 升至 715.5 K，势能降低 22.90 keV；末态没有丢原子，最短 Ti–C 为 1.813 Å。比 v12i 略有改善，仍未通过冷态稳定性，**不启动加热、冷却或凝固**。Fmax 500 eV/Å 仅是紧急异常线。

# Round7 POD/SNAP v12f 拟合与已完成验证归档

本归档建立在分支 `codex/round7-v11-si-hotspot-labels-20261118` 的提交
`beb4c8794c3edbb754d0efc8c00090184ca852b4` 上。该父分支已经保存 4 个 Si 热点标签及锁定划分：
00/02/03 为训练候选，01 为独立开发验证留出。这里另存 v12f POD/SNAP 拟合产物、完整训练/验证构型、相关标签源文件和已完成的筛查结果。

## 内容

- POD v12f 描述文件、系数、FitPOD 输入/日志及误差明细。
- SNAP v12f 模型、系数、参数、FitSNAP 输入/日志及误差表。
- 拟合实际使用的 141 个训练构型（128 个基础构型 + 13 个局部力构型）和 12 个基础验证构型。
- 局部 DFT 标签及审查清单，包括 8 个原始热点标签、2 个 POD 崩塌补充标签；Si 00/01/02/03 已在父分支归档。
- 力列顺序审查、拟合数据清单、Si 01 留出点评估和全模型静态/短时冷态/FIRE 日志及摘要。
- 全模型起始结构没有重复提交：SHA-256 为 `4170c70647fe98c4b1c25960dbae2a557a74adbf48a8a93daf686c8fd16722e7`，与父分支 v10 归档中的 `start_wrapped.data` 一致。

## 目前结果与限制

FitPOD 和 FitSNAP 日志、模型系数文件均已生成。原 `local_training_manifest.json` 中的 “refit pending” 字样已过时；以本归档的系数文件和拟合日志为实际产物记录。

基础 12 构型验证集：
- POD：能量 RMSE 0.0755 eV/atom，力 RMSE 0.6133 eV/Å。
- SNAP：能量 RMSE 0.0533 eV/atom，力 RMSE 0.6521 eV/Å。

Si 01 留出构型的中心力验证未通过：
- POD 向量误差 13.16 eV/Å，方向余弦 0.364。
- SNAP 向量误差 48.01 eV/Å，方向余弦 -0.163。

347,464 原子静态 `run 0`：
- POD Fmax 67.86 eV/Å，F99 9.41 eV/Å。
- SNAP Fmax 269.98 eV/Å，F99 51.37 eV/Å。
- 两模型在该静态点都没有原子超过 500 eV/Å；这不构成长时 MD 或热源准入证明。

已完成的无激光冷态短跑为 1 fs × 10 步（0.01 ps）：固定层 z=0–20 Å，z=20–40 Å 设 300 K Langevin 热沉，其余可动原子按 NVE 积分。末步可动区温度约 2054 K，说明初始结构的势能释放很强；此结果不通过冷态稳定性复核。

原始 FIRE 在第 60 步报告丢失 2 个原子；边界填充副本报告丢失 1 个原子，但该副本还改变了边界、MPI 分区和固定组，不能作严格对照。后续完全按原设置复跑（仅将 lost 策略改为 ignore、限制到 62 步并逐步输出 xyz）时，60–62 步均保留全部 347,464 个唯一 ID，x/y 仅有周期边界坐标轻微超界，z 没有越界。该复跑仍因最大迭代数停止，未收敛；原始 lost-atoms 报错没有被复现为物理逸出，根因仍未查明。诊断输入、日志和汇总已归档；634 MB 逐步 xyz dump 保留在本机，归档记录了其 SHA-256。

SI01 的常规 12 构型总体力 RMSE 仍为 POD 0.613、SNAP 0.652 eV/Å；压缩 Si 热点的独立 01 验证仍局部失效。新 Si ID 244669 构型仅作为训练候选，PBE/PW350/Gamma 标签待收敛和审查，01 保持留出。当前没有通过冷态稳定性验收，也没有启动加热或冷却。

## 校验

归档内 `SHA256SUMS.txt` 覆盖本目录全部归档文件（不包含校验文件本身）。执行：

```sh
cd archive/round7-v12f-pod-snap-artifacts-20261118
sha256sum -c SHA256SUMS.txt
```

## 修正高度模型（单独构建）

另归档了符合 0–200 Å 材料层、200–300 Å 真空区的新几何，458,248 原子，Ag:Ti₃SiC₂ 质量比约 1:1，两相均连续。初始 POD v12f `run 0` 得到 Fmax 22.35 eV/Å。FIRE 和冷态稳定性检查仍待完成；这不是加热准入结论。详见 `models/target_200x200x300_qsgs_material200/`。


## v12i 更新

`v12i/` 保存新增 Ti–C DFT 标签、POD/SNAP v12i 完整拟合候选、Si01 留出评估、两模型全模型静态筛查、POD 冷态 10 fs 日志和压缩坐标，以及复现脚本。详细判定见 `v12i/ROUND7_STATUS_V12I_20261118.md`。


## v12j 更新

`v12j/` 包含完整 v12j POD/SNAP 拟合候选、Si01 留出、新 pair4 DFT 标签、完整模型静态与冷态日志、压缩末态坐标及复现脚本。详细结果见 `v12j/ROUND7_STATUS_V12J_20261118.md`。


## Latest validation: v12k

The current Round7 exploratory POD/SNAP candidate and its validation package are in `v12k/`. See `ROUND7_CURRENT_STATUS_20261118.md` for the current decision and measurements. The package includes the 128/12 data split, 20 local training frames, fitted POD/SNAP coefficients, pair5 DFT label, Si01 holdout evaluation, full-model static summaries, and a compressed 10 fs POD cold-run trajectory. POD cold run reaches about 707 K in 10 fs. The SNAP cold run reaches about 3975 K and creates 2352 Ti-C contacts under 1.8 A; neither model is approved for heating, and SNAP shows severe collapse.


## Latest validation: v12l

The latest POD/SNAP model is v12l. See `ROUND7_CURRENT_STATUS_20261118.md` and `v12l/README.md`. Both models meet the small validation RMSE target, but Si01 remains a difficult holdout and full-model no-laser cold tests reveal instability: POD develops a local Ag-C-C collapse at cumulative 37 fs with Fmax 1399.8 eV/A; SNAP reaches 3330 K and has 1989 Ti-C pairs below 1.8 A at 10 fs. No heating, cooling, or solidification run has started.
