# 纯 Ag 液态 DFT 标签候选

## 构型来源

- 32 个 Ag 原子，三维周期立方盒，边长 8.763314 Å。
- Williams–Mishin–Hamilton 2006 Ag EAM 只用于取样：300 K 平衡、2000 K 熔化、1500 K 密度平衡及 NVT 采样。
- 选取第 44000 步；名义温度 1500 K，LAMMPS 瞬时温度 1384.283 K。
- 采样记录在第 40000 步出现一次 8002 K 的瞬时温度尖峰；选用的 44000 步位于其后 2 ps，EAM 记录温度为 1384 K。它仍是探索性液态样构型，不宣称已证明平衡液体分布。
- 源轨迹以 `dump_modify sort id` 按原子 ID 排序。输入第 i 行映射到源 atom ID i+1；该映射由轨迹末帧坐标核对。
- `candidate_audit.json` 显示密度 8.517 g/cm³、32/32 原子 CNA 为 unknown，且 1 ps 位移分布明显非固态。它们只支持“液态样候选”筛选，不能代替 DFT/MD 相判定。

## GPAW 标签协议

- PBE，平面波截断 500 eV，3×3×3 k 点。
- Fermi–Dirac 展宽 0.10 eV；能量、密度和本征态收敛阈值均为 1e-5。
- GPAW 26.7.0、ASE 3.29.0、四 MPI rank，k-point 并行 4、domain 1、band 1。
- 固定几何单点。输出本征值外推总能与展宽自由能；力按 GPAW/ASE 有限展宽力保存。

## 最终状态（2026-11-23）

GPAW 已于第 20 次 SCF 收敛，运行 10,452.68 s。独立脚本 `verify_label.py` 报告 `PASS_LABEL_INTEGRITY_TRAIN_SPLIT_PENDING`；扩展复核也通过，核对了原子数/顺序/坐标、输入与输出 SHA、有限能量和力，以及 summary 与 extxyz 中 native/free 能量的一致性。完整收敛日志、summary、raw/audited extxyz、label audit、extended review 和每运行文件 SHA 清单位于 `runs/attempt_001/`。大 GPAW 检查点（约 2.36 GiB）未上传，只保留 `checkpoint.json` 元数据。

标签：32 原子 Ag；native extrapolated energy −81.6323852432 eV/cell；force-consistent free energy −81.7484897629 eV/cell；有限展宽力 Fmax 2.84448 eV/Å、F99 2.81687 eV/Å、向量 RMS 1.45580 eV/Å。能量键分别保存，没有混用。

**角色仍为未分配训练候选**，未合并到 `train.extxyz` 或 `dev.extxyz`。结构是单帧 EAM 生成的液态样构型，瞬时温度 1384.28 K、密度 8.5170 g/cm³；文献密度式在该温度给出约 9.1784 g/cm³，该帧约低 7.2%。因此它不代表平衡液态集合；截断能和 k 点收敛也未验证。只能在数据家族复审后决定是否用于探索性拟合，不能用作独立验证或物性准入证据。

`source/` 保存取样输入/轨迹、筛选、GPAW 预检记录和实际标签脚本。软件环境与大 GPAW 检查点不纳入仓库。
