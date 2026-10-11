# D 节点 Ag64 暖机离线热化与结构初筛

日期：2026-10-11 UTC
状态：`OFFLINE_REVIEW_COMPLETE`；不包含新的 MD、DFT 或训练。
任务：`C-D-AGL-1250K-rho097-seed01`

## 结论

D 的 step 0–1000 暖机在归档和恢复层面合格，但热平衡尚未建立。结果归类为：`INITIAL_TRANSIENT = YES`、`THERMALIZATION_NOT_ESTABLISHED = YES`、`NUMERICAL_INSTABILITY = NO_EVIDENCE`、`THERMALIZATION_CANDIDATE = NOT_YET`。现有 0.5 ps、11 个结构帧不能证明液态形成，也不是 DFT 标签。离线审计以 C 基线提交 `50dcc15f6fce508e8d2fbaea4bf63ed7f65b881e` 为准。

推荐：当前保持 HOLD。若用户后续批准，可只评估同一 D 任务从全局 step 1000 延续至 step 2000 的一个 1000 步分段，作为热化趋势诊断；在新的执行器门禁、配置 SHA 和 C 签署完成前，不得执行。现有 `--stage continue` 的终点是 step 120000，不能直接用于短段。

## 来源、归档和状态校验

本审查从 D 分支最新可读取归档独立 fetch、解包并重新核对，不改动 D 分支或原 tar。归档在 D 分支 `cloud-d-b1-dft`，任务归档提交 `11803909226538ab83847d1f06827a8292501139`；随后 D 的恢复回执提交为 `145ee8050cf95a3a1988e868c6b935e10538f723`。

归档对象位于 `codex_warmup_archive/C-D-AGL-1250K-rho097-seed01/b62a818e58a2cf42/`：

| 对象 | SHA-256 / 状态 |
|---|---|
| 原始 tar.gz，100700 bytes | `b62a818e58a2cf423e4cb1b9ed0468ae73ad578a75c4b2cb9fb3f28baf75e7b1` |
| 原始 `archive-manifest.json` | `ab9f0dffa8ef145eac7792dfe1e14e07c1354df1c2152dafec97bad1dc3a57b1` |
| fresh-fetch `REMOTE_READBACK_RECEIPT.json` | `8bc2b17738f5e214098964c0635ef3505134726f8b4d8e7323cd4c5b03591bd4` |
| manifest 内成员 | 18/18 路径、字节数、SHA-256 通过；回执记录 outer/inner SHA 均通过 |
| `thermal.csv` | `aece67f7faddb8f38f379c2ea34967f0ac9b86628d5376d549e40e0e48482140`，6048 bytes |
| `trajectory.traj` | `7b5185a80526cd95d51c07d31f98f3627aa5f4f7242b9f513157e398e68b2818`，56549 bytes |
| `sampler_state.json` | `e08574678b3f8231c0a52ab25febb13c9d06e815b28149eb306d4b451c83ef6f` |
| `source_manifest.json` | `9a694ea8f93655efa7a1c8695669ab2f17520c2d1b90f6df28115b2b1f9c2ad6` |

51 行热数据记录覆盖 step 0–1000（每 20 步一条）；轨迹有 11 帧（每 100 步一帧）。64 个 Ag、周期立方盒边长 10.826824658074 Å、三方向 PBC、任务 seed 125009701 和父族 `PF-AGL-1250K-rho097-seed01` 与 source manifest 一致。原始 D 记录标记 `liquid_QA=NOT_PERFORMED`、`STRUCTURE_PASS=false`、`DFT_LABELS=ABSENT`，数据角色为 `liquid_pilot_candidate_not_training_or_independent_test`。

step-1000 restart 的 `state_arrays.npz` 与末帧轨迹的坐标、动量、晶胞、原子序号、质量和 PBC 逐项完全相等；restart 的 `state_atoms.traj` 坐标和动量也完全相等。checkpoint 保留 `numpy.Generator(PCG64)` 的 packing、velocity、thermostat 三套 RNG 状态，Langevin 版本 5。状态为 `WARMUP_COMPLETE_C_REVIEW_REQUIRED`，不是科学验收状态。

使用的锁定身份：执行器 SHA `6b54d094c26bb9d4a5e4eb5997a31955eb8829a2975042e23fe10527cb943583`；原 C 暖机配置 SHA `2e122b258a48fe54bf603709b5505a301c1257f398c5625eecc2ca94998c0052`；科学参数 SHA `8f792d6206c1b7bae7b900e7ea81445d6cd3e837de2b5979c88bb977650ab7b2`；Ag 势 SHA `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`。势文件本体未在原归档内，本报告也不包含。

归档压缩包本身没有 README.md 或标准 `SHA256SUMS.txt`；原始 manifest 已逐个列出全部 18 个成员的大小和 SHA。此目录的 README、SHA 清单、分析图表和输入副本均是事后补充内容，不是原 tar 成员。`archive-manifest.reference.json` 和 `remote-readback-receipt.reference.json` 是原 D sidecar 的逐字节副本；`inputs/from_D_archive_members/` 是归档成员的提取副本。原 tar、原 D commit 和历史报告没有更改。

## 热日志统计

温度目标列 `target_temperature_K` 在本次暖机为 1600 K；这不是最终 1250 K 平衡段。窗口边界按 global step 含端点，热日志采样间隔为 20 步，即 10 fs。表中标准差为样本标准差；漂移为该窗口对时间做最小二乘直线拟合后的斜率，统一按每 100 fs 表示。

| 时间窗 | 记录数；步数范围 | T 均值 ± SD (K) | T 最大 (K) | T 漂移 (K/100 fs) | 势能均值 (eV)；漂移 (eV/100 fs) | 压力均值 (eV/Å³)；漂移 (eV/Å³/100 fs) |
|---|---:|---:|---:|---:|---:|---:|
| 全程 | 51；0–1000 | 2213.35 ± 650.11 | 4266.30 | −286.17 | −150.094；−5.619 | 0.04645；−0.01792 |
| 前 100 步（含 step 0 和 100） | 6；0–100 | 3229.75 ± 1092.44 | 4266.30 | +5467.44 | −126.817；−73.647 | 0.11478；−0.13913 |
| 最后 200 步（step 800–1000） | 11；800–1000 | 1838.35 ± 105.48 | 2082.21 | −171.67 | −157.922；−2.418 | 0.02061；−0.01295 |
| step 100–500 | 21；100–500 | 2391.48 ± 658.20 | 4062.35 | −897.92 | −148.757；−6.791 | 0.05351；−0.02545 |
| step 500–1000 | 26；500–1000 | 1899.02 ± 144.09 | 2295.43 | −115.20 | −156.414；−1.650 | 0.02568；−0.00572 |

温度在 step 80 / 40 fs 达到 4266.30 K，是日志峰值；step 1000 / 500 fs 的瞬时温度为 1855.25 K，比当时 1600 K 目标高 15.95%。末 200 步均值比目标高 14.90%，且温度、势能与压力的窗口趋势仍在下降。末窗的瞬时波动不大，不能抵消这些均值和趋势仍未稳定的证据。

`sampler_energy_eV` 明确是 `atoms.get_potential_energy()`，即势能，不包含动能。记录器虽读取动能，但没有将它写入此列。`thermal_derived_quantities.csv` 中的动能是按 ASE 温度定义推导值：FixCom 后自由度 `3×64−3=189`，`K=189 kB T/2`；总能仅为“势能+推导动能”，并非日志原始字段，Langevin NVT 下也不要求守恒。step 0、温度峰 step 80、step 1000 的（T，PE，推导 K，PE+K）分别为：

| step | T (K) | PE (eV) | 推导 K (eV) | 推导 PE+K (eV) |
|---:|---:|---:|---:|---:|
| 0 | 1600.00 | −107.864 | 13.029 | −94.835 |
| 80 | 4266.30 | −139.534 | 34.742 | −104.792 |
| 1000 | 1855.25 | −158.721 | 15.108 | −143.613 |

从初始步到 step 80，势能下降 31.67 eV、推导动能增加 21.71 eV，同时温度过冲；到 step 1000，势能累计下降 50.86 eV，压力从 0.14614 降至 0.01581 eV/Å³（约 23.41 至 2.53 GPa）。这支持“随机初态的结构松弛释放势能，造成启动热化瞬态”的解释。它不是严格的单因果证明：轨迹每 100 步才存几何，温度峰 step 80 没有对应结构快照，且 Langevin 热浴同时交换能量。

## Thermostat、积分器和有限尺寸核验

锁定 runner（ASE 3.29.0）在初态先施加 `FixCom()`，再以 `MaxwellBoltzmannDistribution(..., temperature_K=1600, force_temp=True)` 初始化速度；动力学使用 `Langevin(timestep=0.5 fs, temperature_K=目标温度, friction=0.01/units.fs, fixcm=False, rng=PCG64)`。ASE 的 friction 单位为逆 ASE 时间单位，换算后为 0.01 fs⁻¹，特征时间约 100 fs；`dt × γ=0.005`。ASE 文档明确建议小体系用 `fixcm=False` 并单独 `FixCom`，此实现没有重复的 COM 噪声修正。FixCom 移除 3 个平动自由度；记录的 COM 速度数量级约 10⁻¹⁹ Å/fs，符合 COM 已去除。

对 64 原子、189 个动能自由度，在 1600 K 正则分布下单个独立瞬时温度的理论标准差估计为 `1600×sqrt(2/189)=164.59 K`。最终瞬时偏差约 +255 K，即 1.55 个此估计标准差；小体系自然涨落及时间相关性意味着不能单凭末帧偏差宣判失败。当前主要问题是窗口均值持续偏高并继续漂移，而不是末帧单点偏差。只有 0.5 ps 且热浴相关，不能把 51 行当作 51 个独立样本。

压力由 runner 计算为 `−trace(stress)/3 + 2K/(3V)`，单位 eV/Å³，包含动能项。当前为固定体积 NVT 暖机；压力趋势只作结构和热化诊断，不据此自行改盒尺寸或密度。势截断 5.995011 Å 大于半盒长 5.413412 Å；runner 的 ASE 邻居列表支持周期多镜像。这里的 RDF 和局域几何筛查使用最小镜像并将 RDF 限定到 L/2，避免把多镜像成对距离误当作不同最近邻。

## 11 帧几何、RDF 和局域有序度初筛

所有 11 帧均为 64 个 Ag，晶胞和 PBC 一致，坐标有限。按最小镜像统计，704 个“每帧、每原子的最近邻距离”最小值 2.20069 Å，中位数 2.53926 Å，90 分位 2.67680 Å，最大 2.96663 Å；逐帧全局最小距离从 step 0 的 2.20069 Å 到 step 1000 的 2.33547 Å，均高于 runner 的 1.8 Å 几何停止阈值。近邻距离逐渐放松，未见碰撞崩坏迹象。

用全部 11 帧、0.05 Å bin 的纯 Ag RDF，在首壳层识别到平滑第一峰约 2.725 Å（平滑 g≈2.48），第一谷约 3.675 Å。以该谷值作为粗筛 cutoff，平均配位数由 step 0 的 9.53 变到 step 1000 的 11.00；平均局域 Steinhardt Q6 从 0.341 变到 0.367。RDF 峰、配位和 Q6 有变化，但统计来源只有 11 个相互相关的短轨迹帧，且 Q6 cutoff 是从同一小样本 RDF 取出的描述量；它们既不能证明液态，也不能证明晶化。早晚窗口各仅 2 和 3 帧，曲线见图只供诊断。

逐帧几何与逐原子局域量见 `frame_geometry.csv`、`per_atom_local_order.csv`、`nearest_neighbor_distribution.csv`；RDF 数据见 `rdf.csv`。完整 51 行原始热日志逐字节副本为 `thermal.csv`。

## 分类结论

| 分类 | 本次判定 | 依据 |
|---|---|---|
| `INITIAL_TRANSIENT` | YES | 前 40 fs 温度从 1600 K 升到 4266 K，势能和压力快速下降；随机拒绝装填未经势能最小化 |
| `THERMALIZATION_NOT_ESTABLISHED` | YES | step 800–1000 均温仍为 1838 K，窗口 T/PE/P 斜率仍为负，0.5 ps 太短 |
| `NUMERICAL_INSTABILITY` | NO EVIDENCE | 有限数值、轨迹完整、step 1000 checkpoint 与轨迹末态精确一致、距离门未触发、运行正常结束 |
| `THERMALIZATION_CANDIDATE` | NOT YET | 需要更长且更稳定的分段数据；本次不升级任何数据角色 |

“无数值不稳定证据”不等于势函数已通过结构或液态适用性验收。

## 最小追加验证与严格分段门禁（仅设计，当前 HOLD）

当前 step 1000 仍在 `hot_steps=20000` 阶段，所以任何 1000→2000 段仍沿用冻结的 1600 K 热段、0.5 fs、摩擦 0.01 fs⁻¹、原 seed、原盒和原势，不会到 1250 K 冷却/平衡段。checkpoint 中科学参数 SHA 必须仍为 `8f792…ab7b2`，势 SHA 必须仍为 `958f…58447`；不重新随机装填，不重抽 Maxwell 速度，不重置 Langevin RNG。

若以后获批，建议先开发并测试一个明确的显式短段入口（例如 `--stop-at-step 2000`），并让新签署配置绑定任务 ID、owner D、当前归档 tar SHA、step-1000 sampler_state SHA、source manifest SHA、原科学参数 SHA、执行器 SHA 和由 C 在配置文件之外发布的配置 SHA。启动器须把实际配置文件 SHA 与外部固定值比较。要求：

1. 入口必须验证 checkpoint 恰为 step 1000、task/owner 完全匹配、全套状态数组与 PCG64 状态齐全且哈希一致；任何缺失即在计算力调用前 HOLD。
2. 用户/C 新签配置必须显式授权终点 2000 和最多新增 1000 步。CLI 终点与签署配置不一致、没有终点、或试图达到 120000 步一律拒绝。积分循环每步强制检查终点，唯一允许范围为 1000→2000。
3. 原 checkpoint 和 D 原归档只读保留；成功输出写入新 segment 根及新的不可覆盖归档。不得覆盖旧轨迹/日志；step1000 只作连接端点，不重复写入新的积分步。
4. 继续使用原 0.5 fs、1600 K 热段、0.01 fs⁻¹、原势/密度/盒/seed；不得触及 1250 K 目标段、后续119000步或其他五任务。
5. 只允许 1000 步；预估约 200 秒，建议硬墙钟上限 600 秒、本地暂存余量 2 MB（均为计划预算，运行前用实时资源复核）。NaN/Inf、最小距离低于原1.8 Å、资源阈值、阶段超时、checkpoint/归档/远端 fresh-fetch 任一失败立即停止；备份恢复通过之前不允许续段。
6. step2000 后至少对 step1600–1800 与1800–2000 两个非重叠200步时窗比较 T、PE、P 的均值/趋势，并检查几何完整性及 RDF/Q6 的早晚差异；若仍有漂移，继续 HOLD。即使趋势趋稳，此段也只能产生 `THERMALIZATION_CANDIDATE` 供 C 再审，不能声称液态或自动延长。

已附纯 Python `segment_gate_proposal.py` 原型和 `test_segment_gate_proposal.py`。13 项合成单元测试覆盖未签配置、缺失/错误终点、step120000、错误 owner/任务、错 checkpoint SHA、错归档/配置 SHA、错误起点、缺 RNG、科学参数或 source manifest 身份改变、超步数及 step2000 后的逐步停止；合成的 step1000→2000 有效样例仅返回 `VALIDATOR_PASS_ONLY_NOT_COMPUTE_AUTHORIZATION`。测试不调用 ASE calculator、不运行积分、不修改原执行器。当前 `PROPOSED_SEGMENT_1000_2000.json` 状态为 `PROPOSED_NOT_SIGNED`。任何生产实现仍需重新审查代码、冻结新的代码与配置 SHA 并由 C 签署。

## 交付清单

- `README.md`：归档与事后补充文件关系、复现方法和 SHA 说明。
- `SHA256SUMS.txt`：除清单自身外本交付目录全部文件的 SHA-256。
- `analyze_archive.py`、`thermal_windows.csv`、`thermal_derived_quantities.csv`、`analysis_summary.json`：只读分析程序和数值结果。
- `thermal_diagnostics.png`、`structure_prescreen.png`：热日志与结构初筛图。
- `THERMALIZATION_QA_RECOMMENDATIONS.md`：后续热化 QA 判据与停止建议。
- `rdf.csv`、`frame_geometry.csv`、`per_atom_local_order.csv`、`nearest_neighbor_distribution.csv`：机器可读初筛。
- `PROPOSED_SEGMENT_1000_2000.json`、`segment_gate_proposal.py`、`test_segment_gate_proposal.py`：未签署分段方案和无科学计算门禁测试。
- `inputs/from_D_archive_members/`：D 原归档解包成员的逐字节参考副本；不含 Ag 势文件。

本次无 DFT、MD、训练、LAMMPS 或 TTM。D 的其余119000步、其他五条暖机、所有 DFT 和 MACE 训练继续 HOLD；E/F/G/H 不受本报告影响。
