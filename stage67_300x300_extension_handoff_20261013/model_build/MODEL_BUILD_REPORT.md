# Stage67 300×300 Å 横向模型构建报告

## 状态

**BUILD_STATUS = BLOCKED_MISSING_SOURCE_MODEL**。最新指令已明确 `Lx=150.1 Å` 和新热源中心 `(50,150,150) Å`，消除了上一轮的尺寸与中心歧义；剩余阻塞是原始 Stage63 基元与建模脚本缺失。本次没有生成 `300x300_model.data`、LAMMPS restart、可运行 input 或 15×30×30 TTM 映射。没有运行 MD，也没有改势函数、热源功率或液态判据。

## 已核验的原始流程来源

指定 handoff：`cloud-d-stage67-300x300-extension-handoff-20261013`，commit `f4faca36e3987ddc577c49c6ce8c27a1731073da`。其 41 个 ZIP 分卷重组 SHA256 为 `6cbef8dbf96ee2f94bf5c2929b0b622902db12898f879cde00e6c97c56a59547`；包内清单校验通过。

`reports/Stage67_inputs_manifest.csv` 将 Stage67 `model` 输入指向 `Stage63_model_size_test/ttm_run_reproduced_equilibrium/stage63_equilibrium_100.data`，SHA256=`e193a7cde41a06fa04fd32ea9a80c60c4d5089d2a036a0c7fe6957ca57b4813f`。该源文件不在 handoff 包、本地工作区或当前已连接的 GitHub 分支树中。包内只有 200 Å 的 step1650 dump/restart 和带绝对路径的 driver 参考；未含原始 Stage63 data、原始建模脚本或完整可执行输入。清单还引用缺失的 `in.regions_x_surface.inc`。

原始参数可从 Stage67 参考驱动确认：盒子 `150.1×200×200 Å`、PBC `p p p`、timestep `0.1 fs`、TTM `15×20×20`；Ag EAM、TSCBOP、Ag-Ti/Si/C `calibration-v0` 交叉项。势文件只作为历史身份核对，本次未改动。用户指定的目标为 `150.1×300×300 Å`、热源中心 `(50,150,150) Å`，热源功率保持旧设定。

## 可验证的几何数量级

step1650 source dump SHA256=`91d04a6fa3cb1b2fd9eb668b392cad04ecb593620723beac3e12ea4e7b8e85fc`，303,372 atoms，`150.1×200×200 Å`，PBC `p p p`。原子数密度为 `0.0505283145 atoms/Å³`，按原子类型质量映射得到约 `4.71459 g/cm³`。

如果原始 100 Å Stage63 结构通过哈希和几何核验，且原流程确认为从该 100×100 Å 横截面做周期复制，则目标几何应沿 Y、Z 各复制 3 次（LAMMPS `replicate 1 3 3`）：体积与原 200 Å 模型之比为 2.25，计数估算为 682,587 atoms。按 step1650 类型数比例缩放，估算 Ti=238,383、Si=74,646、C=156,654、Ag=212,904。以上是从 200 Å 计数推算，不是已构造或已核实的 target 数量。

该推算不能代替缺失的原始 100 Å 坐标。用 step1650 dump 的某个象限替代 Stage63 平衡结构会把演化后坐标冒充初始模型，因此没有这样做。目标最小原子距离也尚无可测量的 target 坐标。

## TTM 映射

用户指定的 `15×30×30` 网格对应 13,500 cells，尺寸为约 `10.0067×10×10 Å`。旧 `15×20×20` 的 `cell_material_fraction`、`G_mix`、初始 Te 和热源空间文件不适用于目标盒，未复制改名。缺少最终 target data 时无法正确计算格点材料分数、混合耦合、初始 Te 和新热源分布，也无法验证中心 `(50,150,150) Å` 与源功率积分一致。

## 阻塞及恢复条件

需要提供或提交下列原始输入：

1. `stage63_equilibrium_100.data`（SHA256 应为 `e193a7cde41a06fa04fd32ea9a80c60c4d5089d2a036a0c7fe6957ca57b4813f`）；
2. 若原始模型不是直接周期复制，提供 Stage63→Stage67 的实际建模脚本/命令和映射说明；
3. Stage67 manifest 引用的 `in.regions_x_surface.inc`，以及若 restart 必须含完整 TTM 状态，提供同版本 TTM restart/初始化接口。

收到并验证上述输入后，才能按原流程产生 `150.1×300×300 Å` 数据，再生成匹配的 restart/input 和 15×30×30 网格并做成分、密度、最小距离、尺寸、功率积分及哈希验收。
