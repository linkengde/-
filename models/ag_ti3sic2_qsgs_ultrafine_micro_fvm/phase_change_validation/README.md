# 单温 Ag 焓法相变验证包

本目录保存单温隐式焓法的实现、输入参数、解析基准、三维短程算例、数据、图件和资源报告。数值方法验证与实际材料预测分开判断：当前代码/离散验证通过，但Ag物性出处尚不完整，真实电弧边界、液态流动和实际触头凝固预测未通过验证。

## 已运行的验证

- 小型A–H代码检查、潜热/能量守恒、双层界面通量以及独立常物性Stefan熔化/凝固基准。
- 10 μm域上的50³纯Ag、双层和seed 1105 QGS局部数值脉冲—冷却短算例；QGS完成了三档时间步对照。
- 100³纯Ag和双层短算例。100³ QGS在60分钟墙钟限制内未完成，已记录为资源超时；没有启动150³相变循环。

三维热源为仅作用于Ag区域的数值局部体热源，随后六面300 K恒温冷却。它不是电弧热输入，也不能代表触头实际工况。所有三维计算均无液态流动、自由液面、蒸发、反冲压力或表面张力。

## 在当前仓库复跑

```bash
cd /workspace/-/models/ag_ti3sic2_qsgs_ultrafine_micro_fvm/phase_change_validation
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR=/tmp/fvm-phase-mplconfig
python run_small_benchmarks.py
```

完整三维驱动会依次运行多个50³和100³算例，复杂100³ QGS每例最多占用60分钟；仅在确需复跑时执行：

```bash
python run_3d_phase_change.py
```

完成仿真后生成图件：

```bash
python render_validation_figures.py
```

检查本验证包内文件：

```bash
sha256sum -c SHA256SUMS.txt
```

校验和清单以本目录为当前路径，列出本目录中除清单本身和Python缓存外的文件。上游输入文件的独立校验和见 `INPUTS_SHA256SUMS.txt`。

## 主要文件

- `MATERIAL_PROPERTY_AUDIT.md`：参数、单位、适用范围与可追溯来源审计。
- `ENTHALPY_MODEL.md`：焓法、糊状区和离散模型说明。
- `STEFAN_VALIDATION.md`：解析熔化/凝固前沿对照。
- `TIME_STEP_CONVERGENCE.md`：QGS时间步收敛及空间收敛状态。
- `PHASE_CHANGE_VALIDATION.md`：相变验收结论。
- `PERFORMANCE_REPORT.md`：实际耗时、内存与规模判断。
- `results/`：可读CSV/JSON及用于复核图件的计算字段。
- `figures/`：三维形貌、温度/液相率、Stefan前沿、时间步和能量图件。
