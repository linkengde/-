# Ag/Ti₃SiC₂ 势函数工作日志与交接记录

**记录范围：** 2026-10-01 至 2026-10-04（北京时间，UTC+8）

**仓库：** `linkengde/-`，分支 `main`

**基准状态：** 本日志整理时远端提交 `69018884c7b01d7ab3f9c23c4a6b0b7087dd449d`（`Evaluate v12 on expanded interface DFT set`）。

## 工作目的

建立并验证能描述 Ag–Ti₃SiC₂ 界面的多元素势函数，为后续候选结构模拟提供依据。候选结构的 Ti₃SiC₂ 骨架外形应保持；Ag 质量分数以约 59.8775 wt% 为目标，允许小幅调整；界面间距和模型尺寸可按需要调整。最终计划是研究加热、TTM 和压力条件下的 Ag 迁移、熔池及表面变化。

在势函数的界面能量和受力通过独立 DFT 检验前，不把它用于长时间 MD、TTM 或物理结论。

## 前置问题与路线调整

此前的 LAMMPS TTM 10 步诊断证明输入可以短程运行，但不能证明势函数可靠。诊断记录显示 step 0–10 最大力上升，且存在约 1.0002 Å 的 Ag–Ti 近距离接触；高力中保守势力与 TTM 随机力未能明确分离。旧的 calibration v0 使用未经独立验证的 Ag–Ti、Ag–Si、Ag–C 交叉 Morse 项。因此停止延长烟测，转向周期性 DFT 标签和 MACE 多体势验证。

这段前置诊断来自先前对话中的报告；该报告没有包含在本交接仓库中。它说明了为什么选择当前路线，但不作为本日志中 v9–v12 数值的来源。

## 按时间记录

### 2026-10-01：建立交接仓库

- 初始化 GitHub 仓库，作为 GPAW、MACE 数据和模型的跨窗口交接位置。
- 目标从未经验证的 pairwise 交叉项转到“DFT 计算界面能量与受力 → 训练多元素 MACE → 独立评估”。

### 2026-10-03 至 2026-10-04：准备新窗口可复用的计算包

- 整理四个固定几何输入：AgSi 2.2/2.8 Å、AgC 2.0/2.4 Å；每个输入含结构、元素、ID、盒子和中心接触信息，并用 SHA256 清单记录。
- 添加 GPAW PW-PBE 计算脚本、检查点恢复脚本、归档程序及环境说明。
- 固定 DFT 方法：GPAW 26.7.0、PBE、平面波 500 eV、Gamma 点、0.1 eV Fermi 展宽；这些是固定结构单点能量/力计算，不是几何优化。
- 将 v9–v11 模型、v11 数据划分、AgSi/AgC 外部参考和 MACE-MP-0b3-medium 基础模型放入仓库。
- 记录的环境验证：Python 3.12、ASE 3.29.0、gpaw-data 1.2.1、MACE 0.3.16、PyTorch 2.5.1 CPU 版；v11 模型推理成功，`mace_run_train --help` 可运行。完整记录见 `environment_validation.md`。

### v9：力训练基线

- v9 是四元素 Ag/C/Si/Ti 的 force-only 基线，模型名称为 `MACE_periodic_v9_forceonly`。
- 可用交接包没有保存 v9 完整训练日志和全部原始训练设置；不要据此补写未记录的超参数。
- 后续统一评估显示，v9 在五个外部固定结构上的能量 MAE 约 265.2 meV/atom；AgC 外部能量误差约 495.5 meV/atom。它适合作比较基线，不能作为已校准的生产势。

### v10：加入能量参考和 Ag–Ti 能量点

- 在 v9 力/界面数据基础上，加入 fcc Ag、hcp Ti、diamond Si、diamond C 四个 PBE 元素参考；给现有 AgTi 2.3871 Å 帧增加 DFT 能量，并加入 AgTi 2.55 Å 周期界面能量/受力点，保留 2.70 Å 测试点。
- 四元素能量组成矩阵达到秩 4/4，使四种元素能量基准可识别；这不等于所有交叉界面势已经物理验证。
- 统一外部评估中，v10 能量 MAE 约 231.8 meV/atom；AgC 外部能量误差约 300.2 meV/atom。力误差没有稳定改善。
- 数据说明见 `periodic_interface_v4/mace_periodic_v11_energyfocus/data/dataset_manifest.json` 和该目录的 `data/README.md`。其中目录名为 v11，但数据清单记录的是 v10 energy/reference-calibration 数据。

### v11：加强能量拟合重点

- 使用 v10 的 23/2/3 train/valid/test 数据划分：23 个训练结构、2 个验证结构、3 个冻结测试结构；23 帧有力标签，其中 19 帧有能量标签。
- 训练目标强调能量，记录权重为 `energy_weight=100`、`forces_weight=1000`。它沿用了 v10 训练集，并没有补入 Ag–Si、Ag–C 的周期 PW-PBE 界面能量标签。
- v11 在五个外部固定结构上的能量 MAE 降至约 127.5 meV/atom，但 AgC 外部能量误差仍约 203.5 meV/atom，外部力向量 RMSE 约 0.0856 eV/Å。

### 2026-10-04：补齐四个 v12 训练标签、构建并训练 v12

- 完成并归档四个周期 PW-PBE 单点：AgC 2.0/2.4 Å、AgSi 2.2/2.8 Å。逐项检查收敛状态、输入哈希、原子数/元素/ID/坐标、有限能量和力，以及日志和输出能量一致性。
- 用四个 PW-PBE 能量/力标签替换训练集对应的 LCAO 力训练帧；保持原 validation/test 和外部 holdout 分离。数据仍为 23/2/3，训练能量组成秩为 4/4。
- MACE v12 完成 80 轮，选用 epoch 76 checkpoint。模型 SHA256 为 `011280f8d999a62c1603848226fd64c9592740d54e28d4d04f76e2035a034628`。
- 发现并解释训练日志的指标差异：MACE 表格显示的 1800.8 meV/atom 能量 RMSE 把 4 个无能量标签的 force-only 帧当作零能量纳入未掩码汇总；只对 19 个实际带能量的训练帧独立审计，RMSE 为 40.99 meV/atom。这是指标汇总/掩码问题，不是 checkpoint 与评估脚本不一致。
- 三个冻结测试结构上，v12 能量 MAE 为 7.31 meV/atom、力向量 RMSE 为 0.0582 eV/Å；v11 对应为 10.64 meV/atom、0.0632 eV/Å。测试结构只有 3 个，因此证据仍有限。

### 2026-10-04：补做界面距离 DFT 检验并复评

- 新增并验证三个固定几何 PW-PBE 标签：AgC 2.30 Å、AgSi 2.60 Å、AgTi 2.60 Å。输入哈希、结构、收敛、能量和受力均通过归档检查。
- 将三个新点与原有冻结/外部结构一起比较 v9–v12。要注意这些距离扫描点来自已有中间距 motif 的几何扰动，彼此相关，不能视为三个完全独立的界面环境。
- 五个外部结构整体上，v12 能量 MAE 为 9.95 meV/atom，相比 v11 的 127.46 meV/atom 显著降低；整体受力向量 RMSE 基本不变：v12 0.08553、v11 0.08557 eV/Å。
- 分界面看，v12 的 Ag–Ti 力检查改善；Ag–Si 的向量力 RMSE 改善，但中心分离力被低估；Ag–C 的整结构力 RMSE 和中心分离力误差恶化。AgC 2.2026 Å/2.30 Å 两个点的 v12 力 RMSE 分别约 0.1166/0.1197 eV/Å，高于 v9–v11；v12 中心分离力也偏离 DFT。
- 因此报告结论是 `overall_v12_validation = UNDETERMINED`，不是 PASS。外部结构少、几何相关，且生产 MD/TTM 的误差接受阈值尚未定义。

## 当前已完成与未完成

**已完成：** 四个 v12 训练用 DFT 标签及归档；v12 数据集；80 轮训练和 epoch 76 模型；训练指标审计；冻结测试与五个外部固定结构评估；三条新增距离 DFT 标签和扩展评估。证据文件列于下方。

**尚未完成：** v13 数据/训练；足够多、彼此独立的 Ag–C/Ag–Si 界面构型验证；候选大模型静态力验证；无热无压短程稳定性；TTM 场和 checkpoint 连续性复核；稳定熔池；配对压力模拟。没有用 v12 跑生产候选长 MD/TTM，也没有据此宣称势函数可用于生产。

## 下个窗口从哪里继续

1. 先拉取 `main` 最新提交并读取本日志及 `results/v12_validation_assessment.json`，不要重跑已经归档的 DFT 或 v12 训练。
2. 按评估报告建立**独立的 v13 数据目录/划分**，把 AgC 2.30、AgSi 2.60、AgTi 2.60 三个已核验 DFT 点加入训练候选；保留 v12 文件、冻结测试集和未用外部点，不覆盖 v12。
3. 继续为 Ag–C、Ag–Si 补充更独立的界面局部环境/距离构型；对中心分离力偏差较大的点优先取 DFT 能量和力，并保留未参与训练的验证点。不得把同一 motif 的微小距离扰动当成完全独立的泛化证明。
4. 训练 v13 后，统一在冻结测试集及所有预留外部结构上比较 v9–v13 的能量、向量力和中心分离力；报告逐结构结果。若数据仍不足或某交互仍退化，结论标记为未定并继续补点。
5. 只有独立界面受力验证支持之后，才开始候选结构的静态审计和短程力学稳定性测试；通过后再验证 TTM，之后建立熔池和同 checkpoint 的压力/零压力配对模拟。每关通过才继续；未通过则停止。生产模型保持骨架外形及约定 Ag 质量分数，所有结构改动在副本中记录。

## 主要证据文件

- 环境：`environment_setup.md`、`environment_validation.md`
- v9–v11 数据与历史：`periodic_interface_v4/mace_periodic_v11_energyfocus/data/README.md`、`dataset_manifest.json`
- 四个主训练 DFT 点：`periodic_interface_v4/pbe_interface_energy_additions_v12/archive_manifest.json` 及 `calculations/{AgC_d2p0,AgC_d2p4,AgSi_d2p2,AgSi_d2p8}/`
- v12 数据与训练：`periodic_interface_v4/mace_periodic_v12_interface_energy/data/dataset_manifest.json`、`train.stdout`、`results/v12_training_frame_audit.json`
- v12 评估：`periodic_interface_v4/mace_periodic_v12_interface_energy/results/v9_v10_v11_v12_candidate_comparison.csv`、`v9_v10_v11_v12_candidate_comparison.json`、`v12_validation_assessment.json`
- 新增距离 DFT 标签：`periodic_interface_v4/pbe_interface_energy_additions_v12/calculations/validation_distance_scans/`

归档的 GPAW 日志、extxyz 和 JSON 可用于结果复核；模型 checkpoint 及数据文件均由仓库清单记录哈希。不要把“数值可运行”“训练完成”或“能量 MAE 下降”等同于生产级势函数已经通过验证。
