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
