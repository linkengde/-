# 七标签DEV来源审查

审查基准：C 固定提交 `26330e37e3efa2499c6bb304bd5de899e3a6a983`。本文件只记录公开源与来源链，不生成结构、不看封存测试、不改变任何数据用途。

## Ti₃SiC₂ 候选源

本地唯一核实到、可追溯的晶体学文件为 COD 9009647：

- CIF：`research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif`
- CIF SHA-256：`6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3`
- 文献：Kisi, E. H.; Crossley, J. A. A.; Myhra, S.; Barsoum, M. W., “Structure and crystal chemistry of Ti3SiC2,” *Journal of Physics and Chemistry of Solids* 59 (1998), 1437–1443, DOI `10.1016/S0022-3697(98)00226-1`。
- CIF 记载：Ti₃SiC₂，空间群 `P 63/m m c` (No.194)，六方晶胞 `a=b=3.0575 Å, c=17.6235 Å`，Z=2；12原子常规晶胞。
- 现有48原子骨架控制结构组成为 Ti₂₄Si₈C₁₆，是该源的 `2×2×1` 后代；baseline、Ti/C ±0.02 Å均属 `Ti3SiC2_bulk_parent_v19`。

**独立性结论：排除为新DEV源。** 它正是已评分骨架家族，且与既有 Ag–Si 结构提案的底层 CIF SHA 完全相同。不能用另一超胞、重写CIF、位移、晶格变换或不同数据集ID重新命名为独立父族。

尚未核实到第二个实验来源/CIF及可验证样品或文献谱系。COD 和 Crossref 的只读检索在当前环境均被网络代理拒绝（CONNECT 403），因此没有证据可冻结另一来源。所需替代源证据：原始文献或样品标识、CIF原文和校验值、化学式/占位/晶胞、独立来源说明；随后须证明其不是COD 9009647的再收录或同一晶体的转换，并按新源重新映射Ti/C位点。

## Ag–Si界面候选源

现有可审查的候选文件：

- `research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/inputs/AgSi_cod9009647_slab_v19_dev_proposal.extxyz`
- SHA-256：`8fe44bee35e17ebfc433e059925a7a85f0dbb4e775090a75ae1883b3eaa6c36c`
- 52原子：Ag₄Ti₂₄Si₈C₁₆；提案表面终止为Si，Ag单层覆盖；原始输入 `PBC=T T T`，后续A计算明确使用 `PBC=T T F`，因此新包必须锁定结构与计算边界口径。
- 输入晶胞 `c=34.2058765 Å`，原子z包络约`8.0003–26.2059 Å`，沿z合计空隙约`16.0003 Å`；这是旧提案的几何审查值，不能据此认定新界面的真空已收敛。
- 来源字段明确写入 COD9009647 / Kisi et al. 1998，且 manifest 的 `source_cif_sha256` 是上面的旧CIF哈希。
- 该几何被登记为同源数值收敛pilot；后续A公开的三点同样出自此 parent 家族并已评分。

**独立性结论：排除为新DEV源。** Si终止、完整四元素和可读原子组成均可核实，但它没有独立的骨架来源，而且相关几何已进入数值pilot与评分；增加真空、改变间距、扩胞或调整Ag registry仍是该来源后代。

没有找到已核实的第二个 Ti₃SiC₂ 晶体来源及其独立 Si终止界面。候选筛选必须取得来源CIF、明确`(0001)`终止/重构、完整Ag–Ti–Si–C组成、原子映射和PBC，并与已公开训练/开发诊断谱系做父族及周期几何近重复审查。未经独立 evaluator 的盲态谱系回执，不得声称与封存测试互斥；C不会查看封存内容。

## 最小七标签冻结草案与前置条件

目标标签槽位固定为同一个新Ti₃SiC₂父族5点（baseline、一个Ti位点z向±0.02 Å、一个C位点z向±0.02 Å）及同一个新Ag–Si父族2点（固定合法registry下`d0−0.10 Å`、`d0+0.10 Å`）。这七槽位目前均为**来源未选/结构未生成/DFT未授权**，不是可执行输入或已存在标签。

新骨架按来源若具有12原子常规晶胞，可暂估48原子`2×2×1`超胞；若晶胞/基矢或占位不同，需在不改变物理父族前提下重新确定超胞与原子数。不可沿用旧4f atom ID，必须按新来源建立对称位点及ID映射。固体PBC三向周期。已核实旧父族的`10×10×4`仅作参考；新根需固定结构比较倒空间密度网格（建议`8×8×4`与`10×10×4`起始对照），检查native energy、force-consistent free energy及力差，收敛后为5个点锁定同一k网格。

Ag–Si草案目标为Si终止、4个Ag/`2×2`表面胞、足够厚度的完整骨架slab；成本估算暂用80–120原子（中值约100）作为计划包络，不是实测来源或冻结化学式。明确`PBC=T T F`、xy偶极修正；起始面内网格`6×6×1`，同一几何以`8×8×1`复核，依据新表面胞间距调整而不得一律Gamma。初始真空建议至少18–20 Å并增加4 Å同几何对照，检查表面偶极/能量/力变化；不能将旧52原子薄层的设置直接移植。

方法候选：GPAW 26.7.0 / ASE 3.29.0 / PAW 1.2.1，PBE，PW 500 eV，展宽0.10 eV，SCF阈值沿用能量/密度/本征态`1e-5`；每帧保存native energy、force-consistent `free_energy`及全原子三维力，单位eV与eV/Å。全包必须使用一致的能量键；本次没有新计算。
