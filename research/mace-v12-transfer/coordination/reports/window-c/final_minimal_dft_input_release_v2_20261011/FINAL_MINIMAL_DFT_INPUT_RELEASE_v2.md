# FINAL_MINIMAL_DFT_INPUT_RELEASE_v2

**状态：`BLOCKED_BY_SOURCE`。** 本次只核验仓库中既有输入结构的来源、几何身份与 SHA-256；没有启动 DFT/MACE 训练，没有修改结构、数据角色或势函数。该文件不是 DFT 启动授权。

## P1：Ti₃SiC₂ 的 Ti/C 位移输入

几何文件可逐字节核验，但它们来自已评分的 COD 9009647 / Kisi 1998 父族，不满足此前方案要求的独立 DEV 来源。五个历史骨架点（baseline 与四个位移）继续保留 `numerical control` 角色。B 的四个位移 DFT 已完成，不应因为这次输入审查重复计算。

父结构：

- 路径：`research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/inputs/Ti3SiC2_baseline_PROPOSAL.extxyz`
- SHA-256：`157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4`
- 48 原子，Ti₂₄Si₈C₁₆；cell：`[6.115,0,0]`, `[-3.0575,5.2957453441,0]`, `[0,0,17.6235]` Å；PBC=`T,T,T`
- 父来源 CIF：`research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif`，SHA-256 `6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3`
- 来源族：`COD_9009647_Kisi_1998_Ti3SiC2_scored_numerical_control_family`

原子映射以 extxyz `local_new_id` 为准；它是本文件内的 ID，不是旧 LAMMPS ID。

| Wyckoff 位点 | 元素 | 文件顺序（0/1 起） | local_new_id | parent 坐标（Å） |
|---|---|---:|---:|---|
| Ti(4f) | Ti | 2 / 3 | 3 | (1.52876529, 0.88261540, 2.38798425) |
| C(4f) | C | 8 / 9 | 9 | (−0.00001529, 1.76525727, 1.27241670) |

四个已存在的端点如下；同一固定晶胞中只移动映射原子，坐标与输入 SHA 已由当前文件核验：

| 结构 | 位点 / ID | 位移 | 端点坐标（Å） | 输入 SHA-256 |
|---|---|---:|---|---|
| `Ti3SiC2_Ti4f_id0003_z_m020A` | Ti(4f) / 3 | −0.02 Å z | (1.52876529, 0.88261540, 2.36798425) | `5529a040f0f35957f826b78754910970741bff741474e321dc1646a5c0c12bbf` |
| `Ti3SiC2_Ti4f_id0003_z_p020A` | Ti(4f) / 3 | +0.02 Å z | (1.52876529, 0.88261540, 2.40798425) | `5829ebcbb9ae331f6cf1dd7dc5abca0787c07667a55778e3c1eea241938b227b` |
| `Ti3SiC2_C4f_id0009_z_m020A` | C(4f) / 9 | −0.02 Å z | (−0.00001529, 1.76525727, 1.25241670) | `19f1a4011620f5c7f245d0a8c893964aa8ebd1eb647f8a9b6f4547a75cac586f` |
| `Ti3SiC2_C4f_id0009_z_p020A` | C(4f) / 9 | +0.02 Å z | (−0.00001529, 1.76525727, 1.29241670) | `6308961f33ad5fae3f58930570a47044a236e062de243821c20577c78955ff72` |

**P1 判定：** 坐标身份已知；独立 DEV 来源未解决，且四个端点已有 DFT 标签。因此不能把它们作为新的 DEV 输入发布或重标。需要真正独立的新父结构，而不是同一 COD 原型的重命名、扩胞或简单位移。

## P2：Ag–Si slab 输入

现有候选：`research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/inputs/AgSi_cod9009647_slab_v19_dev_proposal.extxyz`，SHA-256 `8fe44bee35e17ebfc433e059925a7a85f0dbb4e775090a75ae1883b3eaa6c36c`。结构含 52 原子（C₁₆Ag₄Si₈Ti₂₄），由 COD/Kisi 同一已评分来源族构建，不具备独立 DEV 身份。

| 项目 | 核验结果 | 放行判断 |
|---|---|---|
| Si 终止面 | 候选中接触侧最上层基底为 Si，平均 z=23.6058765 Å；其下最高 Ti 为 21.58798575 Å。Ag 层平均 z=26.2058765 Å。 | 仅确认该 proposal 的切面确为 Si 顶层，不等于弛豫界面确认。 |
| registry | builder 将 Ag 放在 `(i+0.37)a+(j+0.21)b`（i,j=0,1）；脚本 SHA `ff801ff1545246fa2d22765e9642abe26e69658a10e46cef4b43510c2ac04813`。 | 这是合成坐标偏移，无文献坐标、弛豫或能量证据，不是已确认 registry。 |
| `d0` | 可复现的几何定义候选是顶层 Si 平均 z 到 Ag 层平均 z 的法向距离；现有未弛豫 proposal 的差值为 2.600000 Å。 | 2.6 Å 是 builder 手动给定的初始 gap，不是已证明的平衡/弛豫 `d0`。真实 `d0` 未知。 |
| Ag 层移动 | 方案拟将四个 Ag 原子作为刚性层整体沿 z 平移 −0.10/+0.10 Å，基底、横向 registry 和 cell 固定。 | 这是待执行的几何规则；相应端点结构和 SHA 尚不存在。 |
| slab cell | 向量为 `[6.115,0,0]`, `[-3.0575,5.2957453441,0]`, `[0,0,34.2058765]` Å；候选 PBC=`T,T,T`；当前 z 原子范围 8.00026725–26.2058765 Å。 | 必须保持原件；不能静默改 PBC 或 cell。 |
| 真空 | 沿 z 的周期性空区（Ag 最上层到下一周期基底最下层）为 16.00026725 Å。 | 低于此前方案的 ≥18 Å 目标；+4 Å 真空收敛对照未做。 |

**P2 判定：** termination 可从 proposal 确认；registry、物理 `d0` 和可执行端点未获得来源支持，真空条件也未满足。故 P2 保持 `HOLD`。

## 最终决定

`FINAL_MINIMAL_DFT_INPUT_RELEASE_v2 = BLOCKED_BY_SOURCE`。P1 提供了完整的历史控制几何标识，但不是新的独立 DEV 来源，且 DFT 已完成；P2 缺少可验证的独立父族、可信 registry、弛豫 `d0`、端点文件，并且当前真空间隔低于原定目标。要解除 HOLD，需要补齐新的来源证据和结构文件后再出新版本；本文件不授权 DFT。
