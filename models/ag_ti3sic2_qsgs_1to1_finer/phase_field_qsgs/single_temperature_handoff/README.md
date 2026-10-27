# Ag–Ti₃SiC₂ 单温瞬态热传导—焓法交接包

本目录为独立开展 Ag–Ti₃SiC₂ 双连续 QSGS 结构单温热传导、局部熔化与凝固计算准备数据交接。它复用仓库已有几何、物性候选表、焓法原型和数值验收记录；不复制原子势、不重建 QSGS，也不改变任何原有结果。

## 版本与计算范围

- 交接基线：GitHub `main` 的提交 `1334d51a6f3611cc9fd16c619c433099245990af`。
- 单温模型：热源直接进入 Ag 总晶格焓，使用 Ag 总热导率；**不使用**双温模型的电子热容 `Ce`、电子热导率 `ke` 或电子—晶格耦合 `G`。
- 当前程序是三维连续介质有限体积原型，不是 TTM-MD，也没有原子运动、熔体流动、表面形变或晶粒形成。
- 双温计算及其检查点属于独立运行；不得将双温检查点作为单温初态。本交接包不修改、不覆盖双温任务或文件。

## 目录内容

- `INPUT_MANIFEST.json`：输入文件路径、字节数、SHA-256、格式、用途、最近源提交与核验状态。
- `MATERIAL_PROPERTIES.md`：Ag 与 Ti₃SiC₂ 候选热物性、换算、插值和数据来源缺口。
- `MODEL_SETUP.json`：几何、热源、网格、边界、初态、焓法相变和现有功率统计。
- `VALIDATION_AND_RESULTS.md`：既有单温/双温数值验证、局部输运审计、结果限制与未完成验收。
- `RUN_ENVIRONMENT.md`：复现命令、依赖、当前原型入口及资源记录。

主要代码与几何仍保存在上级目录及 `phase_field_qsgs/`，未复制到此目录。清单列出所有交接所依赖的版本。

## 快速开始：仅复现现有单温短试算

在仓库根目录运行：

```bash
cd models/ag_ti3sic2_qsgs_1to1_finer
python -m pip install -r requirements-ttm-reference.txt
python -m unittest -v test_thermophysical_model test_qsgs_enthalpy_solver test_qsgs_source_integration
mkdir -p runs/single_temperature_handoff
python phase_field_qsgs/qsgs_enthalpy_solver.py \
  --spacing 5A --boundary A --heat-ps 0.1 --cool-ps 0.1 \
  --sample-ps 0.01 \
  --output runs/single_temperature_handoff/caseA_5A_heat0p1_cool0p1.csv
```

脚本会同时写同名 `.json` 汇总。该命令是现有 5 Å、单个 QSGS 周期单胞、Case A 的短试算，不是生产级局部熔池计算。当前 CLI 不能设置原生 Ag 连通图、横向平铺数、热沉、显式 `dt` 或精确二值几何源积分；它使用投影单元中心采样热源。因此不同网格的峰值/源功率不等价，不能仅凭此命令给出收敛的熔池结论。

如需检查已存在的 Stefan 基准：

```bash
python phase_field_qsgs/stefan_benchmark.py
```

## 交接结论

现有文件足以启动和验证单温焓法数值原型。论文级局部熔池或凝固时间的定量结论仍受 Ag 物性原始出处、局部导热路径、热源离散与骨架/外部散热边界假设限制。详见 `VALIDATION_AND_RESULTS.md`。
