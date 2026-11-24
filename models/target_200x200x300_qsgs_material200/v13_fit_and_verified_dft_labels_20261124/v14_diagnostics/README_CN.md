# V14 诊断阶段归档说明

- `V14_STAGE1_DIAGNOSTIC_REPORT_CN.md`：V13 纯 Ag FCC 体积失稳、液态 Ag 与复合体系失稳的诊断结论。
- `v14_stage1_pure_ag_fcc/`：零压弛豫、固定物理密度 NVE、压缩弛豫态 NVE 的输入、日志、轨迹和汇总。
- `instability_replay_v13/`：四个 V13 小体系逐步重放的日志、轨迹、指标和哈希。
- `ag_fcc_pod_static/` 与 `ag_fcc_eos_retry/`：V13 LAMMPS 及一致设置的探索性 DFT/PBE Ag FCC 能量—体积比较。四个 DFT SCF 点收敛，但仅作诊断，不作为正式 V14 训练标签。
- `ag_liquid_precollapse/`：来自历史 V13 验证父轨迹的 0.060 ps 候选；该帧已提取，DFT 正在计算。若加入 V14 训练，父轨迹不可再作为 V14 验证。
- `not_for_fit/`：2.70 Å 过度压缩试算的未收敛日志；不能用作训练或比较。

`run_v14_*` 和 `evaluate_v13_*` 是复现脚本。默认按归档目录相对路径读取 `fit_v13_weight100/`；可用环境变量 `LAMMPS_EXECUTABLE` 指定 LAMMPS 路径。
