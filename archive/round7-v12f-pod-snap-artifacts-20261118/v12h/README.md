# Round7 v12h 候选模型与冷态验证

本目录保存 POD/SNAP v12h 的完整拟合输入、模型、常规验证结果、Si01 留出点结果，以及 458,248 原子模型的静态力和无热源冷态短跑记录。

## 数据和模型

- `candidate/dataset/` 保存 128 个基线训练构型、17 个局部训练构型和 12 个常规验证构型。
- `candidate/POD/`、`candidate/SNAP/` 保存可直接供 LAMMPS 使用的势文件、系数、参数、拟合日志与误差分析。
- `candidate/holdout_si01/summary.json` 保存中心原子 ID 104359 的独立开发留出结果。Si01 始终排除在训练之外。
- `candidate/tic_center_validation/summary.json` 给出三个新增 Ti–C 标签上的中心力对比。这三个局部样本属于训练内检查，不是独立盲测。
- `full_model/` 保存 458,248 原子完整模型的 POD/SNAP 静态摘要和 POD 冷态 10 fs 日志。
- `legacy_small_14007/` 保存旧 Round7 14,007 原子小模型的短跑摘要与日志；它不是新的 200 × 200 × 300 Å 双连续模型。
- SNAP LAMMPS 参数中 `chemflag=1`、`bnormflag=0`，运行时会出现“probably not what you intended”的警告。当前 SNAP 结果按导出的这组参数计算；在继续把 SNAP 当作候选前，需要核对这一设置是否符合本次化学 SNAP 拟合定义。

## 结论

v12h 常规验证力 RMSE：POD 0.604 eV/Å，SNAP 0.685 eV/Å；两者低于探索性全局目标 3 eV/Å。Si01 中心点误差仍是 POD 7.632、SNAP 49.805 eV/Å，热点转移仍未解决。

POD v12h 的全模型静态 Fmax/F99 为 18.016/5.489 eV/Å；SNAP 为 106.049/23.304 eV/Å。完整模型 POD 无热源冷态 10 fs 后无丢原子，但移动区温度从 300 K 升至 814 K，且出现 7 对 Ti–C 间距小于 1.8 Å。它仍未通过冷稳定性筛查。旧小模型的 POD/SNAP 10 fs 温度分别升至约 3,603 K 和 16,895 K。

Fmax 500 eV/Å 仅作为应急上限。通过该上限并不代表弛豫或 MD 稳定。当前不启动加热、冷却或凝固。

最新进度在归档根目录的 `ROUND7_CURRENT_STATUS_20261118.md`。某些脚本使用 `/workspace/potential-fit/round7_pod_snap/` 的绝对路径；新环境恢复时需按该路径放置已归档的完整模型、势文件和数据。
