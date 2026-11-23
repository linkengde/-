# Round7 Ag/Ti3SiC2 界面分析图

本目录包含三类基于已有 Round7 五帧轨迹（0.60–1.60 ps）的界面与 MAX 骨架分析。PNG 便于预览，PDF/TIFF 用于排版；CSV、JSON 和脚本保留数值与方法。

## 图11：界面接触与骨架迁移

- [PNG](Figure11_Round7_interface_contact_migration.png) · [PDF](Figure11_Round7_interface_contact_migration.pdf) · [300 dpi TIFF](Figure11_Round7_interface_contact_migration.tif)
- [时序数据](Figure11_round7_interface_metrics.csv)

固定跟踪 1,787 个热源核心 Ag。核心 Ag 到 Ti/Si/C 骨架的最近距离中位数从 0.60 ps 的 2.83 Å 降到 1.60 ps 的 2.61 Å；距骨架 3.5 Å 内的比例从 68.3% 升至 74.1%。每个核心 Ag 的 Ag–Ti、Ag–Si、Ag–C 近邻数也上升。Ti/Si/C 骨架原子位移按元素和首帧界面邻域/内部拆分。

## 图12：热源核心 Ag 的局部部分 RDF

- [PNG](Figure12_Round7_interfacial_partial_RDF.png) · [PDF](Figure12_Round7_interfacial_partial_RDF.pdf) · [300 dpi TIFF](Figure12_Round7_interfacial_partial_RDF.tif)
- [径向数据](Figure12_round7_interfacial_partial_RDF.csv)

比较 Ag–Ag、Ag–Ti、Ag–Si、Ag–C 四类距离分布。由于结构是成分强烈非均匀的双连续形貌，图中是热源核心 Ag 条件下的局部部分 RDF；长程基线不按均匀体相解释。

## 图13：Ti3SiC2 MAX 骨架局部结构

- [PNG](Figure13_Round7_MAX_skeleton_local_order.png) · [PDF](Figure13_Round7_MAX_skeleton_local_order.pdf) · [300 dpi TIFF](Figure13_Round7_MAX_skeleton_local_order.tif)
- [配位、键长代理及保留率数据](Figure13_round7_max_local_order.csv) · [机器可读报告](round7_interface_report.json)

不套用 FCC/PTM。分别跟踪 Ti–C 和 Ti–Si 首近邻配位、0.60 ps 首近邻原子对的距离变化及其留在固定首壳范围内的比例，并用几何连通分量监测骨架连接。

## 判读与限制

- 该 QSGS 结构为弯曲双连续界面，没有唯一的全局界面法向；因此图11使用局部距离和元素接触，而不是任意选一个 x/y/z 方向画平面浓度剖面。
- 3.5 Å 界面接触、Ti–C/Ti–Si 首壳截断和 3.25 Å 骨架连通都是几何指标，不等同于化学键能或力学强度。
- 1.60 ps 时几何连通最大分量仍约占骨架原子数的 99.8%。这不能证明骨架承载稳定。
- 输入轨迹只有五帧，未保存逐原子力/应力；不能从这些文件判断分叉处应力集中。Round7 开发 MTP 对熔化与界面行为尚未独立验证，运行摘要还记录了 133.006 eV/Å 的最大检查点力。所有图仅描述该探索性轨迹，不证明界面扩散、化学反应、熔化或熔焊。

## 数据和复现

输入为 `codex/upload-round7-fiveframe-trajectory-1p6ps-20261109` 分支固定提交 `878328378869ba4e5860a4081420fcdc59966c52` 中的五帧全体系 dump，以及更正源区分析中的 `round7_1p6ps_source_core_ag.csv`。两份输入没有重复收录在本目录。运行 `make_round7_interface_analysis.py` 需要 OVITO Python、NumPy、SciPy 和 Matplotlib；脚本默认读取 `/workspace/round7-source/...` 中已解包的数据。

计算细节见 [方法说明](round7_interface_methods.md)。
