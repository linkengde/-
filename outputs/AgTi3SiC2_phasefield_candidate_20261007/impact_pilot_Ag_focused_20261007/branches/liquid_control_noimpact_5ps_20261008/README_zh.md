# 100 ps 末态后的匹配无入射对照

从与 `liquid_diagnostic_5ps_20261008` 相同的第446发末态复制起始结构，沿用相同势、边界、固定/热沉分组、300 K Langevin 热沉和 5 fs 高频输出。按25个0.2 ps窗口运行，总时长5 ps；每个窗口重新读入前一窗口的数据并更换 Langevin 随机种子，以匹配轰击分支的分段流程。

本对照不创建新 Ag 入射原子，原子数应保持4,751。用于比较局部温度、Ag 原子运动和局部结构指标在无后续入射时的演变。它与轰击分支的差别是是否继续注入 Ag；所有末态都保存在本目录，不覆盖100 ps或轰击诊断分支。

驱动：`run_noimpact_control.py`；共享分析程序：`analyze_liquid_diagnostic.py`（运行时设 `LIQUID_DUMP_PREFIX=control`）。

这 5 ps 中，坑区 Ag 温度中位数约 2216 K（逐帧统计的 10–90 百分位约 1312–3195 K）；227 个同一批 Ag 示踪原子的漂移校正 MSD 在 1 ps、5 ps 分别为 5.44、31.14 Å²。局部平均 Q₆ 从 0.3875 到 0.3830。与继续轰击支路的数值对比和限制见 [`../liquid_diagnostic_5ps_20261008/REPORT_liquid_diagnostic_5ps_中文.md`](../liquid_diagnostic_5ps_20261008/REPORT_liquid_diagnostic_5ps_中文.md)。
