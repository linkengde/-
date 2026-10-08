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

- 使用交接包中标为 v10 的 23/2/3 train/valid/test 数据划分：23 个训练结构、2 个验证结构、3 个冻结测试结构；逐帧核对确认 23 帧有力标签，其中 **15 帧**有能量标签。19 帧是替换四个 PW-PBE 标签之后的 v12 数量。
- 模型名为 `energyfocus`，但交接包未保留 v11 原始训练命令/日志，不能确认其确切权重。`energy_weight=100`、`forces_weight=1000` 在 v12/v13 脚本中有直接证据。v11 没有补入四个 Ag–Si、Ag–C PW-PBE 界面能量标签。
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

## 上一窗口交接快照（793ec26，v13 开展之前）

**已完成：** 四个 v12 训练用 DFT 标签及归档；v12 数据集；80 轮训练和 epoch 76 模型；训练指标审计；冻结测试与五个外部固定结构评估；三条新增距离 DFT 标签和扩展评估。证据文件列于下方。

**尚未完成：** v13 数据/训练；足够多、彼此独立的 Ag–C/Ag–Si 界面构型验证；候选大模型静态力验证；无热无压短程稳定性；TTM 场和 checkpoint 连续性复核；稳定熔池；配对压力模拟。没有用 v12 跑生产候选长 MD/TTM，也没有据此宣称势函数可用于生产。

## 上一窗口的继续入口（最新状态见下方续记）

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

## 2026-10-04：本窗口续记——v13 已完成训练，新的 DFT 检查仍在计算

- 已拉取并检查上一窗口提交 `793ec26`。原清单 97 个文件的大小和 SHA256 全部匹配。逐帧核查纠正了 v11 能量标签数：v11/v12/v13 分别为 15/19/22。
- **不要重复 v13 构建和训练。** 独立 v13 目录保留 v12 全部训练帧并增加三份距离标签，train/valid/test 为 26/2/3，能量/受力训练标签为 22/26，元素能量组成秩 4/4；valid/test 与 v12 字节一致。
- v13 CPU 训练实际完成 80 轮，选用 epoch 72。模型 SHA256：`85c16e3bf920418001377cf7049a788e201ac50811a4b2b639f40507b8bd5ae8`。日志、选中 checkpoint、数据和评估位于 `periodic_interface_v4/mace_periodic_v13_interface_energy/`。
- v13 冻结测试能量 MAE 7.25 meV/atom、力向量 RMSE 0.05716 eV/Å；一个 C-gap 帧比 v12 变差。两个相关中间距结构的力 RMSE 改善到 0.02412 eV/Å，但与新增训练母结构相关，整体仍 **UNDETERMINED**。
- `pbe_interface_v13_holdouts/` 的 AgC 排布/局部扰动 DFT 已启动四个 MPI 进程并启用 ScaLAPACK；AgSi 随后由同一队列启动。第二窗口不要重复这两点。一次早期未启用 ScaLAPACK 的尝试已保存为中断记录，不是收敛标签。
- 新 AgC 输入最短 1.782 Å 接触实际是 **Ti–C 骨架原子对**，不能误称为 Ag–C 过近接触。既有数据是少原子簇放在约 28 Å 周期超胞中，不能等同于周期延展体相界面；液态/受热/压力环境覆盖仍需补充。
- 用户接受按新数据继续 v14，必要时 v15/v16；暂定每类界面能量 MAE ≤10 meV/atom、力向量 RMSE ≤0.05 eV/Å、标记分离力绝对误差 ≤0.10 eV/Å。数值门槛是初步筛选，样本不足仍为未定。
- 六个新增输入已完成几何、去标签和不重复 v13 训练帧检查。A 在当前两点之后负责两份 v14 AgC/AgSi 盲测；B 负责三份采集点及一份 v14 AgTi 盲测。两个机器各承担四份 DFT，标签逐点回传。
- **开工入口：`coordination/README.md`、`coordination/tasks/*.json`；新窗口先读 `coordination/START_WINDOW_B.md` 并成功领取任务。** 普通 Git push 防止同时领取；各窗口只写分配目录，主日志和总清单由 A 汇总。
- 当前没有生产长 MD/TTM。下一版训练前先归档 v13 新检查结果；v14 预留盲测不进入 v14 训练。v15/v16 根据前一版误差补数据，并保留新的未训练盲测。
- 续记更新（2026-10-04 16:08 UTC）：当前 A 环境 ID 与远端登记一致。AgC v13 DFT 已在四进程下于第 62 轮收敛（能量 −258.955519609678 eV/cell，耗时 5774 s），仍待归档核验；同一队列已启动 AgSi，日志目前到第 1 轮。两份运行输出和检查点留在 `/tmp`，没有重启作业。B 已由另一独立环境 `c5035b48-f43b-4da0-b8e4-e2862817f86a` 领取；16:08 UTC 的远端状态仍显示四个标签待完成。
- 已准备 `mace_periodic_v14_interface_energy/` 的 v13 acquisition 对比、v14 数据构建、训练和盲测评估脚本。它们以核验通过的五个采集点为训练增量、保留三种界面的 v14 盲测，并在评分文件缺失或标签/哈希未通过时停止；尚未运行训练，也未提前标记 v14 通过。
- 为处理 A/B 不同时完成或临时推送失败，A 的 v13 队列现在在两点归档核验后调用发布脚本；脚本只推送已验证的小型结果、更新工作日志/总哈希清单和任务状态，并拒绝检查点及大文件。B 队列原本逐标签自动推送。失败时由原 owner 在原机器恢复普通 Git 同步；另一窗口只消费已经推到 `main` 的核验结果。
- A/B 队列增加独立进度 watcher，每 10 轮将任务名、状态、轮数、耗时和时间戳写入各自的 `coordination/progress/window-*.json`，并同步到任务登记；`sync_tasks.py status --remote` 合并显示两边最近进度。当前 A 作业按旧启动脚本启动，代码推送后会单独接入 watcher，不重启 DFT。
- 续记更新（2026-10-04 16:20 UTC）：A 的 AgC 仍为第 62 轮收敛、待归档核验；AgSi 在本机日志到第 4 轮，远端 heartbeat 最新第 2 轮，四进程作业仍运行。B 的 `AgTi_registry_probe_v13_01` 在 38 轮收敛，但归档距离元数据检查报差值 `2.4538869e-9 Å`，旧阈值 `1e-10 Å`；B 窗口确认结果已保留、计算队列已暂停且当前无计算进程。核对输入后确认差值来自 manifest 在 extxyz 写盘舍入前记录距离；本地六份 v13/v14 输入的最大同类偏差为 `5.22e-9 Å`。三个归档验证器现将仅对比预写盘 manifest 的距离容差调整为 `1e-8 Å`；输入哈希、读回结构/原子位置、能量和 summary 距离检查仍保持严格。B 拉取修正后可对已保留 AgTi 标签重跑归档核验，再从该点继续队列，无需重算。该修复已本地验证输入距离边界，尚待推送及 B 复核结果。
- 续记更新（2026-10-04 16:22 UTC）：B 随后重跑归档核验，AgTi 标签 PASS 并已推到 `main`；远端任务记录显示余下 AgC 标签已启动，另两项仍排队。B 的归档器将 manifest 距离容差设为 `2e-8 Å`（覆盖 8 位小数坐标序列化的最坏界限），摘要距离仍为 `1e-10 Å`，输入哈希与坐标校验未放宽。A 的 v13 和 v14 主窗口归档器同步采用该界限，避免同类输入元数据舍入误差阻塞结果；尚未据此标记任何 A 标签归档通过。
- 续记更新（2026-10-04 16:24 UTC）：从 `main` 独立复核 B 已发布的 AgTi 标签，归档验证记录和四个文件哈希全部匹配；输入哈希、原子 ID/元素/坐标、晶胞、PBC、有限能量受力、摘要能量和 summary 距离均通过。该 DFT 能量为 −176.579846274888 eV/cell，38 轮收敛。A 的 AgSi 本地日志到第 7 轮，仍在运行；A 的两份 v13 标签尚未完成归档核验。

### v13 independent holdout archive

Both assigned lateral-registry v13 PW-PBE holdouts passed the archive checks and were published independently of window B.
- `AgC_lateral_registry_holdout_v13`: SCF converged in 62 iterations; energy -258.9555196097 eV/cell.
- `AgSi_lateral_registry_holdout_v13`: SCF converged in 46 iterations; energy -207.5444919019 eV/cell.
- Evidence: `periodic_interface_v4/pbe_interface_v13_holdouts/archive_manifest.json` and each label's `verification.json`.
- These labels remain excluded from v13 training; score them before considering them for v14.


### 2026-10-04 20:48 UTC: v13 blind scoring and v14 holdouts

- Independently rechecked all four B archive records against their input/output hashes, verification records, summaries, geometries, IDs, finite labels and GPAW logs. All four passed; `AgTi_registry_holdout_v14_01` remains reserved and is excluded from v14 training. No `.gpw` checkpoint is in the archive.
- Ran `mace_periodic_v14_interface_energy/evaluate_v13_acquisition.py` on the five verified acquisition structures plus the frozen v13 test set. v13 is `FAIL_SCREENING`: AgC force-vector RMSE 0.11368 eV/Å and max separating-force error 0.35750 eV/Å; AgSi 0.09037 and 0.16178; AgTi 0.06259 and 0.04392. Energy errors pass the provisional 10 meV/atom screen for all three interfaces, but force thresholds (0.05 eV/Å) fail for all; AgC and AgSi also miss the provisional separating-force threshold (0.10 eV/Å).
- The five new geometries are registry/local perturbations of prior small-cell motifs. This screening result does not establish transfer to extended interfaces, thermal disorder, liquid Ag or pressure. Keep v13 out of production MD/TTM.
- Evaluation artifacts: `mace_periodic_v13_interface_energy/results/v13_independent_holdout_comparison.{json,csv}`. The JSON records v9-v13 model hashes and per-structure errors.
- Started A's two reserved v14 blind DFT labels with four MPI ranks after confirming there were no existing runs or outputs: `AgC_registry_holdout_v14_01` is running and `AgSi_registry_holdout_v14_01` is queued. Their labels stay out of v14 training and model selection. Progress is published to `coordination/progress/window-a.json`.
- Next: finish and verify both reserved labels, build v14 from the five scored acquisition labels, keep all three reserved labels outside training, verify the dataset gates, then train/evaluate v14 against those holdouts.


### 2026-10-04 23:45 UTC: A v14 blind labels complete; B task rebalanced

- A's `AgSi_registry_holdout_v14_01` converged in 47 PW-PBE SCF iterations with four MPI ranks; energy -204.7159976569 eV/cell. The archive verifier passed every check (input/output hashes, geometry, IDs, PBC, finite energy/forces, summary and log); `state.gpw` is not in the archive. Together with the earlier AgC holdout and B's AgTi holdout, all three v14 blind labels are now available.
- Rebalanced the two cloud windows after the original four B DFT labels completed: A owns the remaining DFT result publication, v14 dataset construction/audit, main work log and global manifest. B now owns an independent per-atom localization of v13 force errors and, once A publishes the validated dataset, the single v14 training and blind evaluation run on B's separate compute. Do not rerun DFT or v12/v13 training.
- Added a B residual-localization entry point that ranks atom force residuals, direction mismatch and nearby atoms, and separates errors near the marked interface pair from the rest of each small-cell motif. This is diagnostic evidence for follow-up sampling, not another validation pass.
- The B training entry point will refuse to run unless the v14 31/2/3 dataset, rank 4/4 energy composition, hashes and all three isolated holdout archives pass. A still must run the dataset builder and audit before sending B its second-stage command.


### 2026-10-04 23:45 UTC: A v14 blind labels complete; B task rebalanced

- A's `AgSi_registry_holdout_v14_01` converged in 47 PW-PBE SCF iterations with four MPI ranks; energy -204.7159976569 eV/cell. The archive verifier passed every check (input/output hashes, geometry, IDs, PBC, finite energy/forces, summary and log); `state.gpw` is not in the archive. Together with the earlier AgC holdout and B's AgTi holdout, all three v14 blind labels are locally available.
- Rebalanced the two cloud windows after the original four B DFT labels completed: A owns the final DFT publication, v14 data construction/audit, main work log and global manifest. B now owns independent per-atom localization of v13 force errors and, after A publishes the gated dataset, the single v14 training and blind evaluation on B's separate compute. Do not rerun DFT or v12/v13 training.
- Added a B residual-localization entry point that ranks atom force residuals, direction mismatch and nearest neighbors, and separates errors near the marked interface pair from the rest of each small-cell motif. This is diagnostic evidence for follow-up sampling, not another validation pass.
- The B training entry point will refuse to run unless the v14 31/2/3 dataset, rank 4/4 energy composition, hashes and all three isolated holdout archives pass. A still must build and audit the dataset before sending B its second-stage command.


### 2026-10-05 UTC: v14 dataset built and independently audited

- Ran `mace_periodic_v14_interface_energy/build_v14_dataset.py` after all three reserved v14 DFT holdouts had archive verification PASS. The v14 dataset contains 31/2/3 train/valid/test frames, 27/2/3 energy labels, 31/2/3 force labels, and energy-composition rank 4/4.
- Independently checked train/valid/test hashes against `dataset_manifest.json`; the v13 validation and frozen test files are byte-identical. Confirmed all five acquisition labels were added to training, and all three v14 holdouts have PASS verification, matching output hashes, distinct geometry and blind/excluded roles. The builder's duplicate checks passed.
- B's task is now ready: run the v13 per-atom residual-localization report, then train v14 once and evaluate it against the frozen test and all three independent holdouts. Training is assigned to B's separate cloud CPU; A will review the metrics and update the main log/manifest.


### 2026-10-05 UTC: B residual report reviewed; v14 training underway

- Reviewed `coordination/reports/window-b/v13_force_localization.md`. The largest interface-local residual is AgC strain acquisition (local force-vector RMSE 0.1606 eV/Å); AgC and AgSi registry holdouts also show elevated contact-separation errors. This is diagnostic evidence for follow-up sampling, not a validation pass.
- B completed the localization report and started the single v14 CPU training run. The latest shared heartbeat observed during this update is epoch 22/80 at 2026-10-05T00:19:19Z; independent evaluation has not started. A has no active DFT or model-training process and will review the published evaluation before any v15 decision.


### 2026-10-05 UTC: A prepared residual-targeted v15 DFT candidates

- Reviewed B’s v13 localization report and prepared one fixed-geometry DFT candidate for each Ag-C, Ag-Si and Ag-Ti interface, targeting local force-residual hotspots. The perturbations preserve atom counts and cell shapes, pass the protected-v14 geometry duplicate checks, and keep minimum pair distances above 1.75 Å.
- The inputs and reproducible generator are in `periodic_interface_v4/pbe_interface_v15_targeted_acquisition/`. These labels are explicitly excluded from v14 train/validation/test and the three reserved holdouts. If used for v15 training, new blind geometries must be reserved.
- A is starting the three four-rank PW-PBE single points while B trains v14 on its separate cloud CPU. The first runner attempt stopped before SCF because the task publish path lacked the repository prefix; the scope is corrected and no DFT process was started by that attempt. No production MD/TTM is running.


### 2026-10-05 UTC: v14 screening failed; A starts targeted v15 DFT

- B completed the one 80-epoch v14 CPU run (last optimizer epoch 79; selected checkpoint epoch 71) and published the independent evaluation. Status: `SCREEN_FAIL`. The three blind holdout force-vector RMSEs are Ag-C 0.15367, Ag-Si 0.11782 and Ag-Ti 0.14661 eV/Å, all above 0.05. Holdout energy MAEs are below 10 meV/atom; separating-force errors pass for Ag-C (0.06139) and fail for Ag-Si (0.23490) and Ag-Ti (0.19261) eV/Å.
- On the frozen test, `C_gap_2p9` force RMSE is 0.08263 eV/Å and `Si_gap_2p9` energy error is 15.91 meV/atom, both above their gates. The evidence is still narrow: one reserved geometry per interface from related cluster motifs. See `periodic_interface_v4/mace_periodic_v14_interface_energy/results/v14_validation_assessment.json` and the comparison JSON/CSV.
- V14 is not suitable for production long MD/TTM. A is running three new local residual-targeted PW-PBE labels in parallel with B’s completed v14 cycle. They are excluded from v14; any v15 training will use fresh blind holdouts.


### 2026-10-05 UTC: B assigned v14 residual localization in parallel

- To keep both windows active without duplicating A’s DFT or retraining v14, B is assigned post-hoc per-atom residual localization on the three already-scored v14 blind holdouts. The output will identify local/remote force errors, top-error atoms and neighbors for the next DFT sampling decision.
- This is diagnostic analysis only; the v14 screen remains failed. B must wait for A’s next explicit assignment after publishing the report.

### 2026-10-05 UTC: v14 force diagnosis and evidence-led v15 changes

- V14's evaluated checkpoint used forces_weight=1000 and energy_weight=100; its force loss weight was already ten times the energy loss weight. Raising this weight alone is not the first change.
- The three V14 blind holdouts pass energy MAE (Ag-C 3.23, Ag-Si 4.62, Ag-Ti 1.68 meV/atom) but fail force-vector RMSE (0.154, 0.118, 0.147 eV/Å). Ag-Si and Ag-Ti also fail marked-pair separating-force error (0.235 and 0.193 eV/Å); Ag-C passes that component at 0.061 eV/Å.
- Gap-distance-only undercoverage is not sufficient to explain the force errors: Ag-Si's 2.439 Å holdout is near V14 training distances 2.440/2.463 Å; Ag-Ti's 2.564 Å holdout is near 2.588/2.600 Å. Local registry/coordination and wider force response remain likely contributors, but the limited number of holdouts does not isolate a single root cause.
- B's audit found the three pending v15 residual candidates are small perturbations of V14 training frames (same-ID positional RMS 0.029/0.052/0.044 Å for Ag-C/Ag-Si/Ag-Ti). Treat them as candidate training additions, never as independent v15 validation. The repository has no unused distinct-registry holdout inputs.
- AgC_residual_shell_v15_01 converged in 60 iterations, energy -257.906079 eV/cell, archive verification PASS. The existing queue has advanced to AgSi; AgTi remains queued.
- V15 plan: preserve the active DFT queue; use verified residual labels for training; create and hash new controlled registry holdouts before training; sample Ag-C/Ag-Ti contact shells and broader Ag-Si framework environments; keep V14 optimizer/loss weights fixed in the first data-driven comparison; evaluate force-vector and separating-force errors independently. Do not begin production MD/TTM before fresh tests pass.
- Detailed decisions and limits are in periodic_interface_v4/V15_FORCE_SEPARATION_PLAN.md; B's source-overlap evidence is in coordination/reports/window-b/v15_holdout_design_audit.md.

### 2026-10-05 UTC: additional V14 reference and split-integrity audit

- A frame-level audit of V14 train/test found an exact geometry duplicate: training frame `AgTi_2p7A_LCAO` and test frame `AgTi_2p70_periodic_PW_PBE_energy_force_holdout` have identical 49 atom IDs, elements and coordinates. The training record is marked force-only LCAO; the test record is PW-PBE. Thus the frozen AgTi test is not an independent geometry test.
- At this identical geometry, the LCAO and PW-PBE force vectors differ by 0.08520 eV/Å RMSE across atoms; the largest atom-vector difference is 0.20132 eV/Å. This is above the provisional 0.05 eV/Å force gate. Because the records also differ in calculation protocol/boundary metadata, this comparison establishes a label inconsistency but does not isolate its physical or numerical cause.
- V14 contains four force-only legacy LCAO training frames (two AgC, two AgTi), consistent with 31 train structures but 27 energy labels. Treat their reference compatibility as unresolved until original settings and matched labels are checked. Keep v14 outputs unchanged; do not count its duplicated AgTi test as independent evidence. The v14 cycle already fails on all three newly scored registry force tests, so this finding strengthens the data-integrity diagnosis without changing its `SCREEN_FAIL` decision.
- Consequence for v15: perform exact/near geometry split audits before training; exclude a test geometry if that geometry has any training label; and replace or exclude legacy LCAO force-only labels unless their source and consistency with the PW-PBE target are verified. The v15 plan was updated with these gates.
- The user expects additional historical calculation data to be uploaded. As of origin/main commit `50804fb`, no new uploaded calculation data was present. Revisit provenance and label consistency when it appears; do not infer unmatched LCAO settings from the current extxyz alone.
- The active AgSi v15 residual DFT remains running at SCF iteration 10 with four MPI ranks; the existing queue is untouched.

### 2026-10-05 UTC: review of uploaded V9-V15 force history

- Fetched `c48ae13` and fast-forwarded the A checkout. The CSV has 33 rows (11 structures x V9/V10/V11); all rows map to a source record and the reported energy/force metrics and available or derived separation-force values match their source tables. Its caveat is correct that V9-V11 scores on the later v14 holdouts are retrospective and historical train membership was not reconstructed.
- The accompanying Markdown and CSV each incorrectly say `AgTi_registry_holdout_v14_01` is an exact duplicate of the V14 training frame `AgTi_2p7A_LCAO`. Direct persistent-ID/element/coordinate comparison disproves that: the registry holdout is derived from the 2.70 Å frozen parent but targets a 2.60 Å contact; nearest V14 training geometry has same-ID positional RMS 0.4153 Å and max displacement 0.5746 Å. It is a distinct geometry from that exact train frame, though still correlated by parent-motif lineage.
- The exact geometry duplication is a separate V14 **frozen-test** issue: `AgTi_2p70_periodic_PW_PBE_energy_force_holdout` is exactly coincident with training frame `AgTi_2p7A_LCAO`. The 0.0852 eV/Å LCAO/PW force-vector difference was measured on this frozen-test geometry, not on `AgTi_registry_holdout_v14_01`. Keep these two records distinct when interpreting scores.
- The history report's interpretation otherwise reinforces the main result: force-vector RMSE is far above 0.05 eV/Å across V9-V14 registry checks, even where energy or one central separating-force component improves. The v14 registry holdout table gives Ag-C/Ag-Si/Ag-Ti force RMSE 0.1537/0.1178/0.1466 eV/Å; no interface passes. v14's Ag-Si/Ag-Ti separation-force errors are 0.2349/0.1926 eV/Å. This supports fixing dataset coverage and reference consistency before another training cycle, rather than relying on version number or loss-weight changes.
- The uploaded report was not edited in B's path. This correction is recorded here and in `V15_FORCE_SEPARATION_PLAN.md` for A/B handoff.
- AgSi residual DFT is still active at SCF iteration 20, four MPI ranks; AgTi remains queued.

### 2026-10-05 UTC: AgSi v15 residual label verified and published

- AgSi_residual_shell_v15_01 converged in 49 PW-PBE SCF iterations on four MPI ranks. The archive verifier passed input hash, atom count, element/ID/coordinate identity, cell, PBC, marked-pair distance, finite energy/forces, energy agreement, and convergence-log checks. An independent ASE read confirmed exact input/output coordinates and IDs and finite 51x3 forces. Energy is -208.5784241289574 eV/cell; the marked Ag-Si distance is 2.4505465 Å and separating force is 1.9719203 eV/Å.
- Compact result, summary, progress, GPAW log, launcher log, partial manifest, task status, and the synchronization fix were pushed to `main` at `9c5c670`. No `.gpw` checkpoint was included. The state is `AgC` and `AgSi` acquisition labels complete; `AgTi_residual_shell_v15_01` remains queued.
- Publication first failed because the task's `allowed_publish_paths` included a future v15 registry-holdout directory that did not yet exist; `git add` rejected the missing pathspec. `coordination/sync_tasks.py` now stages only allowed paths that exist. The queue script also skips labels already marked completed, preventing repeat archive/publication churn on resume.
- Reviewed B's registry proposals: AgSi and AgTi pass its geometry and duplicate screens, but are rigid registry shifts of existing small-cell training motifs. They may test local registry interpolation only; they do not establish transfer to new morphology or thermal states. The AgC proposal has a 1.7413 Å C-Ti pair in the unchanged framework. B's blanket all-atom 1.75 Å cutoff flags it, but this pair needs species- and connectivity-aware comparison before judging it physically unreasonable. A needs a valid Ag-C holdout geometry as well as the AgSi/AgTi proposals before v15 training.
- No v15 dataset construction, training, MD or TTM has started. Next: resume only the queued AgTi DFT; ask B to audit the C-Ti framework pair and search for an Ag-C candidate under species-aware geometry checks; freeze/hash approved blind holdouts before training.

### 2026-10-05 UTC: species-aware Ag-C review and v15 holdouts frozen

- B published its species-aware Ag-C geometry review at commit `ac415b5`; A fetched `main` and independently checked the 57-entry checksum inventory. C#604587–Ti#604590 recurs across eight framework geometries at 1.7372–1.8102 Å with the surrounding C-Ti shell, supporting a locally compressed framework first-shell bond rather than an isolated Ag overlap. A blanket 1.75 Å all-pair cutoff was not appropriate for this pair.
- The Ag-C registry candidate passes B's empirical same-motif species-pair screen and exact/near duplicate audit against 27 same-formula, same-ID references. It is the same geometry as B's earlier unlabeled proposal, translated relative to its v14 training parent; its novelty is local Ag registry only. Its SHA256 is `f528a051246a1a89546da6f55ee3e3b66222cd08caefca039d77989049d25135`.
- A froze AgC, AgSi, and AgTi registry candidates as label-free v15 blind inputs under `periodic_interface_v4/pbe_interface_v15_registry_holdouts/`. Source frame IDs/species, cell/PBC, rigid Ag-only translation, unchanged framework coordinates, and absence of energy/force labels were checked. All eight entries in the input/script SHA256 inventory pass. Inputs remain excluded from v15 training, validation and model selection; their DFT labels are not yet computed. These probes are correlated local-registry tests, not independent morphology coverage.
- AgTi_residual_shell_v15_01 is active with four MPI ranks. Its latest published heartbeat is SCF iteration 10 at 2026-10-05T04:52:35Z; the GPAW log shows iteration 10 and energy -177.041803 eV/cell. It has not converged or been archived. B's geometry-only task is complete; no DFT/MACE/MD/TTM was run by B.
- No v15 dataset, training, model scoring, production MD or TTM has started. Next: allow the active AgTi residual job to finish; archive/verify/publish it; run the three frozen registry holdout single points serially; then audit the v15 train/holdout split and legacy LCAO force-label compatibility before training.

### 2026-10-05 UTC: A reviews B inherited-method tracing

- Synced B completion `36d5db5` and independently verified all 21 entries of its SHA256 list. The report at `coordination/reports/window-b/v15_inherited_method_provenance_trace/report.md` corrects the prior missing-source_method classification: 14 retained rows have alternate `dft_method` declarations; one remains unresolved. None has hash-linked original calculation/convergence/software evidence. Serialized label lineage is confirmed, original-run provenance is not.
- Integration requirements: preserve both method fields with explicit declaration-versus-verification status; never fill absent parameters from config names. Keep energy/free-energy conventions distinct. Different declared k meshes are a convergence-review issue, not proof of a common model-error cause.
- No automatic blanket exclusion: B reports excluding all 15 gives composition rank 3/4. A must review the unresolved AgTi 2.387 row first and document any exclusion/relabeling decision; the 29-frame scenario retaining rank 4/4 does not resolve other provenance gaps. No production dataset or training was authorized by this trace.
- AgC and AgSi V15 blind labels are archived and verified (AgSi published `a3276c4`, 47 SCF iterations). The queue has started AgTi with four active MPI ranks. B tracing is complete; no additional calculation has been started by B. Formal V15 dataset integration and training are pending.

### 2026-10-05 UTC: A integrates limited V15 29-frame screening data

- Synced B revision `8069026`; selected explicit unresolved-AgTi exclusion independently of blind-label values. A build completed with 29/2/3 frames, finite labels and composition rank 4/4. Historical valid/test bytes and blind-input isolation checks passed. The new production screening builder pins B draft SHA256 and refuses nonempty data.
- Training provenance remains limited: 14 partial method declarations, 12 historical explicit declarations not reverified, three verified residual archives. This is screening data, not proof of unified reference-method consistency or model validity.
- AgTi blind DFT is active (latest inspected iteration 6). Training remains pending until the DFT queue completes. B is assigned isolated training/evaluation entry preparation only, with no blind-label access or model execution.

### 2026-10-05 UTC: A reviews B continuous-queue outputs and integrates entries

- Independently checked 24+41+36=101 SHA256 entries across training/evaluation draft, batch-priority and reference-convention publications; all passed. Read each README/report. Six proposed future labels avoid selecting both members of any of the 11 documented near-pairs; geometry selection remains a proposal requiring rescreening before DFT.
- Adapted entry output boundaries to new children of the V15 directory, preserving original B artifacts, all input pins and numerical settings. An initial default preflight rejected A's appended README because it was pinned; restored the original README and moved entry usage to ENTRY_USAGE.md. Default preflight then passed 29/2/3 counts, rank4 and frozen-input isolation without reading blind outputs. Actual training/evaluation not executed; script runtime remains unverified.
- Convention review supports current residual REF energy linkage to GPAW extrapolated/native energy, not free energy. Historical native aliases are lineage, not original-run verification. Preserve REF targets; energy offsets alone do not establish the cause of force errors.
- Latest inspected AgTi blind log is iteration32, unconverged. Training waits for the DFT queue and all archive checks; no production MD/TTM.

### 2026-10-05 UTC: V15 trained and first blind screen FAIL

- Training launcher exit0 confirmed; all epochs0-79 and Done present. MACE selected epoch74 by train/valid criteria. Frozen model SHA256 f258bee17bc46e9399c1f1d13bb5ed30d4a03997b33baac53b5e58c2df6ebc1a; selection record uses one-based75 plus zero-based74 explicitly.
- Reviewed B metric checks and fixed immediate derived-metric finiteness checks and selection-snapshot hash binding before evaluation.
- First fresh registry screen completed: AgC energy0.8571meV/atom, vectorRMSE0.13528eV/A, separation0.38744eV/A; AgSi4.6344,0.31328,0.03630; AgTi1.3162,0.11751,0.35013. All energy gates pass, all force gates fail; only AgSi separation passes. Overall FAIL; no MD/TTM. Different V14/V15 probes are not a matched direct improvement comparison.
- Published training logs/configs/models/checkpoint, frozen selection and evaluation/hash evidence. B assigned post-evaluation force localization with explicit inference authorization; no retraining or DFT by B yet.

### 2026-10-05 UTC: V16 wider-span acquisition owners assigned

- Reviewed B execution pack and V15 force-localization report; hash checks passed. AgSi marked Si error is primarily transverse, so a scalar separating-force pass is insufficient.
- Integrated unchanged six candidate geometries into separate A registry and B distance acquisition directories, three labels per registered environment. Reused proven four-rank runner/archive/queue guards, recording native and free energies separately. A explicitly authorizes these assigned DFT jobs; no V16 training or model claim. B can consume the standing queue without user per-task forwarding.
- This batch samples complementary axes but does not independently displace the full high-error framework or provide new unseen validation. V16 requires future withheld geometry design, label-role ledger and dataset audit before training; do not promote old scored V15 probes to fresh validation.

### 2026-10-05 UTC: B assigned matched force/separation driver diagnosis

- User requested direct investigation of error drivers before subsequent changes. Added priority analytical job after existing B DFT, before validation geometry design. Compare V14/V15 on common already-scored geometries, decompose longitudinal/transverse and framework responses, and rank coverage/reference/optimization hypotheses with discriminating tests. No claim of causal identification from simultaneous dataset changes; no blind V16 outputs, new DFT or retraining in the diagnosis.

### 2026-10-05 UTC: V16 active queue and independent B label check

- Confirmed A registered instance and four active MPI ranks; AgC registry log reached iteration30, not yet converged. Existing queue/checkpoints retained without restart.
- Independently rechecked B AgC distance compact archive: PASS and convergence59 iterations, every declared member/input hash, ordered IDs/elements, positions/cell/PBC, finite energy/force arrays and summary agreement, and converged-log evidence all pass. No method/force-weight changes made.
- B matched-geometry driver diagnosis remains a separate future job, not the already completed V15 localization. New withheld V16 geometry design remains pending; training cannot start merely because six acquisition labels exist.

### 2026-10-06 Beijing: A reviews matched-driver report and B V16 labels

- Synced diagnosis publication and verified all ten package hashes from its declared relative working directory. Same-geometry V14/V15 vectorRMSE: AgC0.142882/0.135282, AgSi0.335893/0.313283, AgTi0.116872/0.117514; all separation errors worsened. These simultaneous dataset edits do not isolate one causal factor.
- Strongest observed error mode is transverse/local-neighbor response, especially AgSi. Provenance uncertainty remains a risk without proof it caused these probe errors. Metric/order/PBC checks on these references passed. No evidence justifies blindly increasing force weight/epochs; preserve controlled recipe for the next reviewed data cycle.
- Independently rechecked all three B V16 distance archives: declared hashes, convergence, geometry/ID order, cell/PBC, finite forces and summary energy pass. A AgSi registry calculation is active (latest inspected SCF5), AgTi queued.
- B is designing new withheld geometries. Added subsequent geometry-only paired transverse/framework training-acquisition design, with roles and proposed scales explicitly separated; no extra DFT or model run authorized under that design job.

### 2026-10-06 Beijing: three new V16 local-registry holdouts frozen for B DFT

- Reviewed B final candidate manifest/contact audit and package hashes. Independently checked identity, unchanged framework, label absence and no near overlap against V15 training, six V16 acquisitions and scored V15 holdout inputs using separate Ag/framework RMS rule.
- Preserved three candidate bytes/hashes in pbe_interface_v16_blind_holdouts, assigned B four-rank DFT. Frozen before V16 model training; all labels excluded from training/validation/model selection. These remain inherited small-cluster motifs, not independent morphology or thermal tests.
- A AgTi final training acquisition remains active; no V16 training yet. B now has executable DFT/archive tasks via standing queue without per-task user forwarding.

### 2026-10-06 Beijing: A AgTi interrupted job recovery

- Actual process inspection found no GPAW ranks; original AgTi log stopped at SCF25 and repository running heartbeat was stale. Reported interruption explicitly and preserved the iteration20 state.gpw (1.87GB), all original logs and input.
- Confirmed same A instance and GPAW26.7 runtime. Launched four-rank recovery into a new `_resume01` directory using saved GPAW parameters, with restored atom numbers/positions/cell/PBC and XC guards. Original outputs were not overwritten or recalculated from initial density. Recovery script/record retained for audit; convergence and archive verification still pending.

### 2026-10-06 Beijing: V16 dataset and controlled training entry ready

- Verified all nine archive hashes/PASS chains. V16 built from retained V15 training29 plus six reviewed acquisition labels, giving35/2/3 and rank4. Labels finite, source/identity checks and new blind-input isolation pass; historical valid/test copied byte-for-byte. New local-registry test labels not read by builder or used for model selection.
- Adapted V15 entry to V16, pinned new inputs/manifest/foundation and kept80 epochs/seed45/lr1e-4/energy100/force1000. Full preflight including archive integrity passes. This controlled cycle tests reviewed data augmentation; provenance/motif restrictions persist, no production MD/TTM.
- B assigned independent V16 integration audit without fresh test-label access; A proceeds to training on idle CPU after DFT completed.

### 2026-10-06 Beijing: V16 metadata repair, completion review and first screen

- B audit found provenance summaries omitted six V16 additions. Fixed builder and current manifest maps to totals35, preflight count guards, entry text and policy-script pins. Preserved before-repair metadata/pins; train/valid/test labels/bytes unchanged. Full preflight PASS.
- All epochs0-79, Done and selected epoch79 verified. Launcher session exit receipt unavailable; recorded null instead of fabricating0 and required explicit full-epoch/Done/checkpoint matching before evaluation. Frozen selected model SHA918aad24341c6892d65ec0d31b19baf610ac0efde1cfd9d2e1c0567320e0cf98.
- First withheld-local-registry screen FAIL: AgC energy0.6979meV/atom, vectorRMSE0.12950eV/A, separation0.04170eV/A; AgSi3.0214,0.09615,0.07928; AgTi1.6120,0.10529,0.32266. All energy gates pass; all force gates fail; AgC/Si separation pass, Ti fails. V15/V16 test geometries differ so these are not matched improvement deltas.
- B assigned metadata repair check and matched V15/V16 error localization on already-scored geometries, with inference explicitly allowed but no tuning/new DFT. No long MD/TTM.

### 2026-10-06 Beijing: targeted paired diagnostic acquisition approved4+4

- Reviewed and independently checked eight label-free geometries from B targeted design: six framework-contact ±0.04A points plus AgSi transverse Ag ±0.15A. Exactly one framework atom changes with Ag fixed, or rigid Ag-only transverse shift; IDs/cell/PBC preserved. Input hashes checked.
- A owns four positive signs and B four negative signs; B first completes matched V15/V16 diagnosis. These near-parent correlated points test controlled local sensitivity, not independent model validation. Same PW-PBE reference recipe/native-free convention retained.
- No automatic incorporation/training: review paired DFT response, conflicts and exact geometry roles first; new V17 unseen validation must be designed separately before training.

### 2026-10-06 Beijing: A accepts evidence-based V17 recipe after cumulative review

- Reviewed B V9–V16 report and hash inventory10/10 PASS. Independently recomputed all8 model aggregates from48 CSV rows; V9/V15/V16 matched pooled vectorRMSE0.2905/0.1852/0.1464 confirm49.6% V9 reduction and21% V15 reduction. Per-interface gates remain failed.
- Adopt unchanged V16 recipe for first V17 data-only cycle; finish8 signed diagnostic labels and review reference sensitivity/identity before integration. Reject blind force-weight/epoch increases and passing scalar separation as complete force validation. Current evidence supports local/transverse coverage as strongest actionable hypothesis, not sole proven cause.
- A decision/evidence arithmetic saved in V17_REVIEW_DECISION.json. New truly withheld geometry-family evidence remains required; parent-correlated paired labels do not provide independent morphology/thermal coverage. A AgC positive framework log latest observed37 SCF; B negative four-point DFT queue started.

### 2026-10-07 UTC: V17 paired DFT progress

- A AgC framework-neighbor positive label converged in58 SCF iterations; archive verifier passed input hash, convergence, four-rank run, atom/element/ID/coordinate identity, cell/PBC, central-pair distance, finite labels, summary-energy agreement and converged-log checks. Published in `2a20057`; the large `state.gpw` checkpoint was not included.
- A AgSi transverse positive label converged in47 SCF iterations. Its compact archive passed the same checks and is published in `2a20057`. Native energy is -208.7862420599 eV/cell, free energy -209.7169782357 eV/cell, marked-pair separating force +1.820831 eV/Å, and maximum atomic force 4.804412 eV/Å. This is a near-parent diagnostic/training acquisition, not an independent validation result.
- Independently checked B's completed negative-sign archive: all four records in `pbe_interface_v17_parallel_targeted_acquisition/archive_manifest.json` have `status=PASS`, matching per-file SHA256 values, `scf_converged=true`, and `verification.json` PASS. B's AgC, AgSi transverse, AgSi framework and AgTi points are complete.
- Thus A has2/4 and B4/4 labels verified. A's AgSi framework and AgTi positive points remain queued. After all eight matched labels are verified, B's assigned paired-response analysis and split-geometry design still precede V17 dataset construction/training. No model training, production MD or TTM has started.


### 2026-10-07 Beijing: V17 paired diagnostics complete; validation roles reviewed

- A's four positive-sign diagnostic labels are archive-verified: AgC framework 58 SCF, AgSi transverse 47, AgSi framework 49, and AgTi framework 38. The AgTi compact archive and task completion were pushed at `43b250f`; its elapsed time was 5804 s. B's four negative-sign archives were independently rechecked against their per-file SHA256 manifests and convergence records; all eight diagnostic labels pass. These remain correlated training diagnostics, not validation.
- Reviewed and accepted B's nine geometry-only proposals: three development-validation points owned by A and six withheld-test points assigned to B. Both published geometry/package SHA256 inventories pass. The candidate audit reports no duplicate/near-duplicate hits against its recorded inventory and role overlap. All candidates still inherit small-cluster parent motifs; morphology and thermal-state independence are not established.
- Created a hash-bound A authorization for only the three development-validation labels. A's static source/input/identity preflight passed. The six withheld-test labels are not authorized and cannot be accessed until a V17 selected checkpoint and selection record are hash-frozen.
- B's paired-response analysis is now assigned and unblocked. A's three four-rank development-validation PW-PBE labels are queued. No V17 training/evaluation has begun; the next training decision follows B's signed-response report.


### 2026-10-07 Beijing: corrected V17 validation runner before SCF

- The first A AgC development-validation launch passed static authorization/geometry preflight, then stopped at import before SCF because GPAW 26.7 requires `gpaw python` for this MPI mode. No DFT label was produced. The launcher traceback is preserved in `development_validation_attempts/AgC_v17_dev_validation_01.launcher_attempt_01.log`.
- Patched the guarded execution-pack runner to use the selected venv Python for ASE preflight and the matching `gpaw python` CLI for MPI calculation. The four-rank GPAW communicator smoke check passed; the pinned PW-PBE settings and all geometries remain unchanged. Re-running static hashes and input preflight before restarting the A queue.


### 2026-10-07 Beijing: finalized GPAW MPI argument forwarding

- The corrected launcher exposed a second pre-SCF CLI issue: GPAW parsed driver options as its own. Preserved this second traceback as `development_validation_attempts/AgC_v17_dev_validation_01.launcher_attempt_02.log`. No DFT label or SCF directory was produced by either failed attempt.
- Added GPAW's `--` argument delimiter. A four-rank `gpaw python` smoke script confirmed all ranks see MPI size4 and receive the driver arguments; the development-label authorization/hash/geometry preflight passed again. Both candidate package SHA256 inventory and A role-decision hashes pass.
- B has started the now-unblocked V17 paired-response analysis on the eight verified diagnostics. A is restarting the three authorized development-validation labels with the corrected runner.


### 2026-10-07 Beijing: A development-validation DFT resumed

- After two logged pre-SCF launcher failures, the corrected GPAW CLI/argument forwarding passed a four-rank smoke check and the exact role/input preflight passed. A's AgC development-validation PW-PBE run is now active on four MPI ranks; local GPAW log reached SCF iteration1. AgSi and AgTi remain queued. No output label exists yet.
- B's task card reports its V17 signed-pair response analysis running. Withheld-test DFT remains blocked pending a frozen V17 checkpoint selection.

### 2026-10-07 UTC: V17 response review and development-validation queue

- Reviewed B's completed eight-pair response report. Its archive checks pass; three V14 parent geometry-hash declarations do not reproduce under the V17 canonical hash, while exact dataset/frame/ID/cell/PBC identity and direct child-displacement checks pass. Native/free-energy secants differ for three pairs, so the V17 target convention must be explicit before dataset construction. Full evidence is in `coordination/reports/window-b/v17_paired_response_analysis/report.md` and `analysis_final/analysis.json`.
- A's three development-validation labels are running serially. AgC is active on four MPI ranks at SCF iteration30 at 2026-10-07 06:25 UTC; no convergence/archive result yet. AgSi and AgTi remain queued. Two pre-SCF launcher failures remain documented and produced no labels. The six withheld-test labels remain sealed.
- Assigned B `v17_energy_convention_and_parent_hash_audit` through `coordination/tasks/window-b.json` (published in commit `274f9c7`). B will trace source evidence and hash algorithms only; it must not read A's dev outputs or any withheld labels, edit production data, or run calculations.
- V17 has not been trained. After all three A labels pass archive checks, integrate only after reviewing both B findings and development-validation metrics; then train/evaluate V17. No production MD/TTM.

### 2026-10-07 UTC: AgC development label archived; final-checkpoint hang corrected

- The AgC development-validation point converged in58 SCF iterations. Its energy/force summary and extxyz were complete and finite; archive verification passed input hash, convergence/log, atom/ID/coordinate/cell/PBC identity, finite labels, method and summary-energy checks. The compact archive is under `pbe_interface_v17_main_targeted_acquisition/development_validation_archives/AgC_v17_dev_validation_01/`; no `state.gpw` was archived.
- After compact outputs were written, the original runner remained in final `calc.write(state.gpw)` for about20 minutes. Source review showed that collective write was called only on rank0. Preserved the files, stopped the stuck ranks after SCF completion, documented the termination in `launcher.log`, and archived without repeating the DFT. Updated the execution pack to omit this unnecessary final checkpoint and record future launcher commands/exit status. The periodic SCF checkpoints remain in the local run directory.
- Refreshed and verified both execution-pack and A role-decision SHA256 inventories. AgSi is next; AgTi remains queued. No V17 training or held-out label read has occurred.
- B's parent audit is hash-verified and supports native extrapolated energies for the three traced V14 parents, with free energies separately retained; legacy parent hash algorithm remains unreproduced. Assigned B `v17_training_pool_provenance_audit` to assess all35 V16 train rows plus the eight published training diagnostics, reusing earlier evidence and keeping A dev/withheld outputs sealed. Task/log updates are published in the current main snapshot.

### 2026-10-07 UTC: B source audit verified; remaining dev labels resumed

- Synced B's `dacb5cb` source-audit report and verified every artifact hash in `v17_energy_convention_and_parent_hash_audit/SHA256SUMS.txt`. The three V14 training parents use native extrapolated `REF_energy`; free energies remain separate. Their legacy frame hashes cannot be reconstructed, while row identity and all eight signed diagnostic displacements pass. This conclusion is limited to those three parents; inherited pool consistency remains unknown.
- AgC development-validation archive verification passed and was published. A's old runner stalled after the compact outputs while rank0-only code called collective final `calc.write`; preserved the outputs and did not repeat SCF. Removed the unnecessary final write and added launcher command/exit logging. Execution-pack and role-decision SHA256 manifests both verify.
- Resumed only the remaining AgSi and AgTi labels with four-rank guarded execution. AgSi passed preflight and is active; its SCF iteration has not yet appeared in the current heartbeat. AgTi remains queued. Six withheld labels remain sealed; no V17 training has started.
- Assigned B `v17_training_pool_provenance_audit` to reconcile the 35 current V16 training frames plus the eight published V17 training acquisitions, reusing old provenance reports and excluding all A dev and test outputs. This will inform V17 inclusion/energy-convention decisions before data construction.

### 2026-10-07 UTC: automatic B queue replenishment

- The remote B task card shows `v17_training_pool_provenance_audit` running at iteration 0 on its registered instance. No second concurrent task was started.
- Added a standing A duty to replenish B's queue before the current task finishes, with one active B job at a time and explicit dependency/release conditions. B continues the next unblocked queued task without requiring the user to forward each assignment; a repository update still cannot wake a stopped session.
- Queued `v17_split_leakage_audit` after the current provenance audit. It checks geometry-role overlap among the 35 V16 train geometries, eight paired V17 training inputs, three development candidates, and six withheld candidates. It is restricted to geometry-only evidence and cannot read held-out labels or outputs.

### 2026-10-07 UTC: clarify B audit startup status

- B reported that `v17_training_pool_provenance_audit` had only been registered: no row ledger or report had been produced. The task card's `running / iteration 0` was registration status, not substantive audit progress.
- Clarified the immediate first milestone in `window-b.json`: publish a hash-checked source inventory and 35-row ledger skeleton, then proceed through the full 35-row/eight-label review. The following geometry-only split audit remains queued.

### 2026-10-07 UTC: reviewed B V17 audit outputs and assigned anchor follow-up

- Verified both B SHA256 manifests: all 10 provenance-audit artifacts and all 4 split-audit artifacts match. Reviewed reports/ledgers. B found 21/35 V16 rows directly archive verified, 14 partial, and all eight paired V17 training diagnostics directly verified; the complete 43-row pool remains not fully original-run-certified. The four unary references are among the unresolved rows; archive-supported-only rank falls to3.
- Split audit covered52 geometries and405 cross-role pairs, with zero exact/near matches under stated geometry thresholds. It reports shared parent-family correlation separately; this does not establish broader morphology or thermal coverage. No blind/development outputs were opened in B's work.
- A AgC and AgSi development-validation labels pass archive verification; AgTi is still running on four MPI ranks, last inspected at SCF iteration26.
- Replenished B's now-empty queue with `v17_element_reference_anchor_recovery_audit`, focused on the four unary reference rows and the smallest evidence-backed way to preserve a rank-4 composition basis. No DFT or training is authorized for B.

### 2026-10-07 UTC: V17 development DFT complete; B rank-recovery analysis assigned

- AgTi V17 development-validation DFT converged in41 SCF iterations on four MPI ranks. Its compact archive SHA256 checks and `verification.json` all pass; native energy is -174.3254015209843 eV/cell and free energy is retained separately. The large state.gpw was not included.
- All three development-validation points now have compact archive verification PASS: AgC58, AgSi47 and AgTi41 SCF iterations. These remain development evidence, not the six sealed withheld tests.
- Reviewed B's unary-reference audit: original records for Ag4/Ti2/Si8/C8 were not found in searched history; exact stored geometry/label lineage is preserved, but method/version/k-points/SCF/getter remain unverified. Removing all four unary rows from the full43-row pool leaves rank4; the 29-row archive-supported core remains rank3. The report correctly limits this to composition rank and makes no quality claim.
- Assigned B `v17_certified_core_rank_recovery_analysis` to enumerate the smallest provisional/relabel subsets that raise the 29-row directly supported core to rank4. This is metadata-only and does not read A development or any withheld outputs.

### 2026-10-07 UTC: frozen V16 development score on V17 structures

- Scored only the three hash-verified V17 development-validation archives with the frozen V16 checkpoint (selected epoch 79; model SHA256 `918aad24341c6892d65ec0d31b19baf610ac0efde1cfd9d2e1c0567320e0cf98`). The existing MACE CPU environment loaded the model and completed inference successfully. No six withheld-test labels or outputs were opened.
- Result is `DEV_SCREEN_FAIL`. Native-energy absolute errors are 2.918 meV/atom (AgC), 2.845 (AgSi), and 1.570 (AgTi), all under the provisional 10 meV/atom threshold. Force-vector RMSE is 0.13764, 0.07289, and 0.10153 eV/Å, all above 0.05; marked-pair separating-force absolute errors are 0.20439, 0.20347, and 0.18565 eV/Å, all above 0.10.
- Largest per-atom force errors occur on Ti ID 604600 (AgC, 0.4590 eV/Å), Si ID 366936 (AgSi, 0.2035), and Ti ID 479027 (AgTi, 0.3154). These three related development motifs are a narrow screen, not a broad independent transferability test.
- Artifacts: `periodic_interface_v4/pbe_interface_v17_main_targeted_acquisition/evaluate_v16_development_validation.py` and `development_validation_evaluation_v16/` (per-structure JSON/CSV, per-atom CSV, SHA256 list). Next: review B's minimal rank-recovery analysis and settle uncertain-frame/energy-anchor policy before building or training V17. No production MD/TTM.
- Follow-up localization of the same frozen-V16 residuals shows the AgC 4.5 Å contact shell RMSE is 0.1632 versus 0.0596 eV/Å outside; AgTi is 0.1253 versus 0.0650. AgSi is 0.0805 versus 0.0591, with high transverse errors on Si#366936, marked Si#366934 and Ag#685172. AgC's largest residual is transverse Ti#604600; AgTi's are marked Ti#479027 and neighboring framework C#480805. This directs the next data check toward paired framework-neighbor and transverse-registry coverage, while keeping the three development labels out of V17 training.
- Published a reproducible script/report/JSON/SHA list under `pbe_interface_v17_main_targeted_acquisition/development_force_error_localization_v16/`; all listed hashes verify. Expanded B's task card to request an immediate rank-recovery heartbeat and an automatic report-only integration-options follow-up. At this log time, B had not yet published a start event; Git cannot wake an ended B session.

### 2026-10-07 UTC: V17 training active and B queue replenished

- The fixed four-thread V17 provisional run `training_provisional_02` is active. Epochs 0–76 completed by 14:41:53 UTC with finite logged train/validation metrics; full 80 epochs remain required by the frozen recipe. No test split or withheld labels were opened.
- Published A progress at epoch 76. B reconciled the four existing archive records at 13:47 UTC without rerunning calculations; as of 14:00 UTC, the approved frame0 DFT still has no start heartbeat; the B task card explicitly advances to it, and one unified resume prompt has been provided. Its queue now specifies, in order: reconcile the four existing B archive records (no DFT), execute the single approved frame0 relabel, independently audit A’s committed V17 training dataset without reading development reference values or withheld outputs, and publish the V9–V17 review after A freezes and scores the model. These repository assignments do not wake an ended B session.
- This V17 dataset remains provisional: 14 inherited rows have partial original-source evidence, and the selected model’s three development scores are not a blind validation. Production MD/TTM remains blocked.

### 2026-10-07 UTC: B V17 work verified; first training export rejected by frozen epoch rule

- Verified B archive reconciliation, frame0 DFT and independent dataset audit. The frame0 result is 32 atoms, C8Ag10Si2Ti12, PW-PBE plane-wave 500 eV/Gamma/Fermi 0.1 eV, 50 SCF iterations, native energy -212.00151021835651 eV and free energy -212.58193242090846 eV; compact SHA256 list and verification PASS. The source row declared 1x4x2 k-points; Gamma convergence is unresolved.
- B audit PASS for the unchanged provisional43 pool: 29 directly archive-supported and14 partial-source frames, exact rank4, finite labels, no exact/near training-to-dev geometry overlap; withheld outputs remained unopened. Its audit report and DFT archive hashes were checked. Frame0 is not yet integrated. B is now preparing an isolated certified30 core-plus-frame0 builder draft; it does not edit A data or read dev reference values.
- `training_provisional_02` completed epochs0-79 and the runner reached its final receipt, but MACE exported the lowest-validation-loss checkpoint at epoch72 because the runner omitted checkpoint-retention flags. The frozen V17 scorer requires selected epoch79. This model has no selection record, was not scored with `evaluate_v17_dev.py`, and is ineligible for use. Its status, logs and artifact hashes are preserved under that run directory.
- Corrected `run_v17_training.sh` to save every checkpoint and fail unless MACE exports epoch79. Updated the builder to require/pin B’s complete aggregate archive manifest. Rebuilt the provisional dataset; train and dev extxyz SHA256 values are unchanged, manifest now pins the complete archive, exact rank4 remains, and preflight passes. Started `training_provisional_03` at 14:53 UTC with the same data and optimization settings; only checkpoint retention and the epoch79 export guard were added. Epochs0-79 completed by15:55:31 UTC with finite metrics; the runner exported the fixed epoch79 model and reached its final success receipt. The 3-frame development-only evaluation then completed with FAIL: energy error 2.755/1.771/1.221 meV/atom; force-vector RMSE 0.1283/0.0613/0.0942 eV/Å; separation error 0.0948/0.1909/0.2135 eV/Å for AgC/AgSi/AgTi. Thus all energy gates pass, all force-vector gates fail, and separation passes AgC only. Model, epoch79 checkpoint, full logs, selection record, metrics, evaluation and SHA256 lists are prepared for publication. No test split or withheld label was opened.


### 2026-10-07 UTC: V17 fixed-epoch screening completed; B review unblocked

- `training_provisional_03` completed 80 epochs on four CPU threads. The selected model is the predeclared epoch79 artifact; stdout records epochs0–79, `Training complete`, `Done`, and successful runner receipt. Model hash `559b12c024fe8e2155e26ffc7ca48cf9fa1c2682ea37ff9c5367af14a5cc43d3`.
- Frozen development evaluation (`evaluation_dev_01`) is **FAIL** under the unchanged gates. AgC: energy 2.755 meV/atom, force RMSE 0.1283 eV/Å, separation error 0.0948 eV/Å. AgSi: 1.771, 0.0613, 0.1909. AgTi: 1.221, 0.0942, 0.2135. Energy passes all; force fails all; separation passes AgC only. These three previously used dev structures are not fresh blind tests.
- The selection record, selected model, epoch79 checkpoint, full training logs/metrics, development CSV/JSON and SHA256 inventories are saved. V17 remains screening-only; no production MD/TTM.
- B's `v17_iteration_review` is unblocked and prioritized to run first; the report-only certified30 core draft follows. A will use that review before freezing V18. No withheld labels were opened.


### 2026-10-08 UTC: V18 clean30 controlled comparison integrated

- A read and SHA-verified B's V17 review and clean30 proposal (main16bbe99). Accepted the recommended fixed-recipe comparison before changing optimizer weights or adding coverage labels. V17 on the same three dev geometries improves all vector averages but remains FAIL; AgTi separation/max regression remains a target.
- Integrated the proposal byte-for-byte into `mace_periodic_v18_clean_core/data/train.extxyz` and copied the same V17 dev bytes to valid.extxyz.30/3/0 frames, exact rank4, finite labels, source hashes and train/dev geometry isolation pass. All30 rows have direct archive support; coverage removal and Gamma-versus1x4x2 frame0 caveats prevent a sole-cause inference. No withheld labels opened.
- Prepared the same seed45/medium/batch1/lr1e-4/wd5e-7/energy100/force1000/80epoch/final79/four-thread training recipe. The runner guards final epoch and automatically writes a frozen selection plus dev/per-atom diagnostics after successful training. Disk guard protects against incomplete checkpoint storage.
- B is assigned an independent V18 entry audit, then frozen V17 vector diagnosis with explicitly limited already-scored development inference, then V9-V18 review when the new model/results are published. No repeated user task forwarding is required for an active B session.

- V18 training_clean_core_01 launched at 2026-10-08T00:01:40.972560+00:00; preflight and initial disk guard PASS,4 CPU threads. This event is a real process launch; epoch progress will be recorded separately.

- B V18 independent entry audit PASS (affae08); all source/count/rank/role/byte checks pass. A accepted the selection/data snapshot binding recommendation and strengthened the evaluator, and records actual OMP thread count in selection. Launcher status still derives from the intended successful set-e runner flow; standalone finalizer is not an independent exit receipt. Pair Cartesian projection is retained for these existing marked pairs, not generalized to arbitrary periodic pairs. B frozen V17 vector diagnosis started at37b3571.

- A verified B frozen V17 diagnosis at196cbfd/88f7bb0 and accepted the signed-pair explanation: AgTi marked Ti projected force error worsens from-0.203264 to-0.232128 eV/A, while markedAg changes only-0.017617 to-0.018657; this explains separation-error regression despite reduced bulk vector RMSE, without identifying a sole training cause. AgC largest residual remains transverse Ti604600; AgSi markedSi/Ag differential projection remains high. B is now assigned label-free residual-direction candidate preparation during V18 training; new DFT requires A review and a separate owner assignment.

- V18 epochs0–16 completed by00:10:18 UTC with finite metrics and about4 active CPU cores. About2.5GiB remains; remaining checkpoint footprint is covered. Added V18 artifact directory to A explicit publish scope. B geometry design is queued; no new DFT has been launched.

- Extended this V18 finalizer to publish curated completed artifacts and update A status/worklog after fixed79/model/data/evaluation checks. It refuses unrelated local/index changes, verifies owner, uses explicit files (not all80 checkpoints), and stops on ordinary push/rebase failure without reset/force. This makes the current training+scoring+publication stage self-contained; later scientific recipe decisions still require A/B review.

- V18 epochs0–25 completed by00:14:44 UTC; metrics finite,4 CPU threads. Successful completion pipeline now includes curated ordinary-Git publication with owner/clean-index guards; only final epoch79 checkpoint is pushed. Intermediate checkpoint hashes are explicitly marked local-only. B residual-direction geometry design remains ready pending its start heartbeat.

- V18 training log now confirms all80 epochs, Training complete/Done and epoch79 export (00:44:21 UTC). Original tool process/launcher exit receipt is no longer available; no selection/evaluation/summary had been generated, so automatic completion was not confirmed. Recover without retraining: require exact exported-model/epoch79-checkpoint tensor equality and complete logs/data hashes, explicitly record original training_exit_code=null, then score only the same3 development frames. No metric gates or split policy changed.


### V18 completed controlled screening (UTC 2026-10-08T01:11:29.366489+00:00)

- Clean30 fixed-recipe run completed80 epochs; selected/exported epoch79 verified. Frozen same3 development screen: FAIL. No withheld labels opened.
- AgC: energy 5.522354 meV/atom; vector RMSE 0.135904 eV/A; separation error 0.213996 eV/A.
- AgSi: energy 0.096705 meV/atom; vector RMSE 0.057137 eV/A; separation error 0.190539 eV/A.
- AgTi: energy 3.146385 meV/atom; vector RMSE 0.095893 eV/A; separation error 0.228859 eV/A.
- Final model/checkpoint, selection, full logs, artifact hash manifest and per-atom evaluation saved; only selected checkpoint is published. B V18 review dependency satisfied. This remains a reused-development/provisional comparison with coverage and frame0 mesh confounds.

### A review: V19 six residual-direction proposals
Independently verified six recommended candidate hashes, finite label-free geometry, element/ID/pair identity, cell/PBC, exact moved-atom displacement and species-pair minima (coordinate serialization tolerance 3e-8 A). Optional lammps_type is absent; no types invented. Geometry screening PASS does not establish physical transferability or independent validation; parents are already-scored development families. DFT not started. Reproducible review at mace_periodic_v18_clean_core/review_v19_candidates.py and .json. B V18 review is unblocked by published frozen model/results; next queued task is a launch-disabled six-label execution pack, planned A positive/B negative signs (3 each), pending official owner registration.

### A accepts B V9–V18 cumulative review
Reviewed v18_iteration_review/report.md and verified all four SHA256 entries. Matched V17/V18 scores and per-pair contributions agree with published development results; B independently verified all 99 frozen epoch79 state tensors. Coverage removal and Gamma frame0 addition confound attribution: source certainty alone is not a force-error mechanism. V19 will prioritize AgC transverse Ti and AgTi Ti/framework differential response, preserve existing recipe/gates, and prepare independent validation before fitting. B v19_dft_execution_pack is ready now; completed review must not be repeated. Actual DFT remains unstarted pending execution-pack verification and owner registration.

### V19 A input staging and storage blocker
Prepared three positive-sign geometry-reviewed candidates byte-identically in pbe_interface_v19_main_residual_acquisition/inputs with SHA256 manifest. B execution pack is registered running (design-only). No DFT started. Disk has only369MB available; previous runner writes full state.gpw. A identified untracked, unselected V17/V18 epoch0–78 checkpoint files as a concrete cleanup proposal, preserving selected79/models/logs/all tracked artifacts. No files deleted; deletion requires user approval.

### V19 execution integration
User-approved deletion completed for158 untracked nonselected V17/V18 epoch checkpoints (~4.99GiB); selected79/models/logs/tracked files retained. Free disk5.4GiB. Verified B execution-pack SHA256 entries and integrated separate production pbe_interface_v19_residual_acquisition; report drafts remain disabled. Explicit A positive/B negative owner registration, four-rank MPI/ScaLAPACK and actualCPU checks, no-active-run guard, same-state unfinished-run stop, exact geometry/hash/finite/convergence archive verification and compact per-label publication. Training and withheld labels remain untouched.

### V19 active queue follow-up
A first AgC positive acquisition has four live MPI ranks (about100% CPU each); initialization active, no SCF iteration or error yet. B negative queue assigned but latest remote progress still execution-pack completed, so B launch not asserted. Added dependency-gated six-label signed-response analysis then independent V19 dataset-entry audit to B rolling queue. Previously scored V17 development parents are now acquisition families and must not serve as fresh V19 validation; independent validation remains required before fitting/selection claims.
