# V14 诊断阶段归档说明

- `V14_STAGE1_DIAGNOSTIC_REPORT_CN.md`：V13 纯 Ag FCC 体积失稳、液态 Ag 与复合体系失稳的诊断结论。
- `v14_stage1_pure_ag_fcc/`：零压弛豫、固定物理密度 NVE、压缩弛豫态 NVE 的输入、日志、轨迹和汇总。
- `instability_replay_v13/`：四个 V13 小体系逐步重放的日志、轨迹、指标和哈希。
- `ag_fcc_pod_static/` 与 `ag_fcc_eos_retry/`：V13 LAMMPS 及一致设置的探索性 DFT/PBE Ag FCC 能量—体积比较。四个 DFT SCF 点收敛，但仅作诊断，不作为正式 V14 训练标签。
- `ag_liquid_precollapse/`：来自历史 V13 验证父轨迹的 0.060 ps 候选；DFT 已收敛，原子映射已核对，V13/DFT 力差和最短 Ag–Ag 对径向力已记录。它仍是训练候选，不是已加入的 V14 拟合数据；若加入，整个父轨迹不可再作为 V14 验证。
- `legacy_ag_cooling_v13_force_checks/`：冻结 V13 对已有 6 帧 Ag 熔化/冷却 DFT 标签的外部力诊断；历史标签为 PW350/Gamma 且能量经过常数平移，因此未直接并入拟合。
- `ag_u3_reference_force_checks/` 与 `ag_u3_reference_eos/`：Foiles 等人 Ag_u3 EAM 在已有 Ag 轨迹和 FCC 体积扫描上的纯银参考检查。该单元素势适合 Ag 采样/对照，不能单独代表四元复合体系。
- `V14_TARGETED_DEVELOPMENT_PLAN_CN.md`：以 Ag 熔化、液态演化、冷却凝固为目标的分阶段开发与验收范围。
- `not_for_fit/`：2.70 Å 过度压缩试算的未收敛日志；不能用作训练或比较。

`run_v14_*` 和 `evaluate_v13_*` 是复现脚本。默认按归档目录相对路径读取 `fit_v13_weight100/`；可用环境变量 `LAMMPS_EXECUTABLE` 指定 LAMMPS 路径。
