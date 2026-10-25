# 主路线准入判定：Ag–Ti₃SiC₂ 原始交叉势模型

日期：2026-10-25。按“先试主路线，未通过则切换 B”的执行顺序，对仓库原始 QSGS 模型及 `Morse calibration-v0` 做了当前静态复核，并结合此前已归档的冷态短跑结果判定是否进入加热。

## 判定

**主路线当前不通过加热准入，停止在冷态门槛；后续切换到几何约束型 Model B 做探索性开发。** 这不是熔池加热失败的结论：本次没有施加热源，也没有启动熔化模拟。原因是原始模型无热源时仍有明显的力学/温度漂移，而所需的生产级 TTM-MD 耦合尚未实现。

## 本次主路线静态复核

- 输入：`in.layered_static_check`，读取原始 917,378 原子 QSGS 数据和 `potential.interface_calibration_v0.inc`。
- 边界：`p p f`；只运行 LAMMPS `run 0`，无最小化、动力学、恒温或热源。
- 结果：LAMMPS 22 Jul 2025 Update 4 完成，耗时约 5 s；最大原子力为 3.9877933 eV/Å，Newton 层最大力同为 3.9877933 eV/Å。完整日志保存在 `primary_route_static_run0.log`。

该值与仓库此前的静态诊断一致。此前已归档的无热源动力学筛查也没有建立可用冷态起点：例如 K=10 eV/Å² 弱约束 25 fs 后活动区约 339 K、Ag 约 341 K，撤除约束后的 Ag/骨架最大力约 4.277/4.274 eV/Å，骨架最大位移约 0.538 Å；刚性骨架 25 fs 对照中 Ag 温度升至约 331.6 K，Ag 最大力约 3.97 eV/Å。细节见 `cold_pilot/RESULTS_COLD_K10_0P025PS_20261025.md`、`cold_pilot/RESULTS_RIGID_SCAFFOLD_0P025PS_20261025.md` 及 `LAYERED_MODEL_PLAN.md`。

## TTM 启动条件

当前已编译 LAMMPS 包含原生 `ttm/mod`，但它不能直接满足本模型的 `p p f` 自由表面、Ag 相掩码及“横向高斯 × 深度指数”源项要求。仓库里的 `ag_ttm_reference.py` 是数值参考原型，尚未接入 LAMMPS 原子积分；它不是可运行正式轨迹的生产实现。因此本次没有用不匹配的 TTM 设置或直接原子加热来替代目标方法。

## 切换到 Model B 的范围

Model B 已单独归档在 [`geometry_constraint_B/README.md`](geometry_constraint_B/README.md)。它把骨架作为固定几何障碍，以 WCA 纯排斥项替代 Ag–骨架 Morse 项。B 的静态力筛查比主路线改善，但 100 fs 无源 NVT 末温约 184 K、仍未达到 300 K 冷态验收，且有 63 个 Ag 原子不在最大连通分量。因此，切换到 B 表示后续探索路线选择，不代表 B 已获准进行正式 TTM 熔池计算；其势参数也只是探索性障碍参数。

## 复核命令

从 `models/ag_ti3sic2_qsgs_1to1_finer/` 运行：

```bash
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp \
  -in in.layered_static_check \
  -log /tmp/ag_primary_route_static.log \
  -screen /tmp/ag_primary_route_static.screen
```

该命令仅复核静态结构和力，不进行 MD 或加热。
