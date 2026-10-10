# 液态 Ag 密度与候选采样器来源核验（只读研究笔记）

日期：2026-10-11（Asia/Shanghai）；未运行 MD、AIMD、DFT、训练或模型推理。

## 密度结论

**RHO_REF_1250K_VERIFIED = NO。** 1250 K 是约比常用 Ag 熔点 1234.93 K 高 15.07 K 的近熔点目标温度；熔点数值在本次仍未取得 NIST 页面正文，保留为公认表值待原文复核。此前的约 9.3 g/cm³仅为熔点附近量级，不能升级为 1250 K 的审定密度。没有可核查的 rho(T) 方程、表号、测量压力和不确定度，因此六个盒尺寸仍阻断，不能生成实际坐标。

正常公开访问测试：
- https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440224&Mask=4 — tunnel 403；未取得正文。
- https://api.crossref.org/works?query=liquid%20silver%20density%20temperature&rows=5 — tunnel 403。
- https://www.google.com/search?q=liquid+silver+density+temperature+equation+1250 — tunnel 403。
- https://doi.org/10.1103/PhysRevB.33.7983 — tunnel 403（EAM论文，不是密度证据）。
- https://www.ctcms.nist.gov/potentials/entry/1986--Foiles-S-M-Baskes-M-I-Daw-M-S--Ag/1986--Foiles-S-M--Ag--LAMMPS--ipr1.html — tunnel 403。
- https://openkim.org/id/EAM_Dynamo_FoilesBaskesDaw_1986_Ag__MO_947112899505_001 — tunnel 403；此具体 OpenKIM 标识仅为检索尝试，正文和存在性均未核实，不能引用为已确认模型。

下一步来源审查要求：取得 Iida & Guthrie, *The Physical Properties of Liquid Metals* (1988) 或等价同行评审实验文献的银液态密度表/拟合正文；记录版次、页/表号、适用温度/压力、原始单位、方程参数与不确定度，并在 1250 K 内插（不得无说明外推）。该要求属于科学输入准备，不是本次生成/计算授权。

64 Ag 立方周期盒的确定公式：
`V_A3 = 64 * M_Ag / (N_A * f * rho_ref_gcm3) * 1e24`；
`L_A = 22.547972740316116 * (f * rho_ref_gcm3)^(-1/3)`，使用 M_Ag=107.8682 g/mol、N_A=6.02214076e23 mol^-1。
三组 f=0.97,1.00,1.03 的 L/L_ref 分别为 1.01020478645,1,0.99019544705。公式/相对比率可以先准备；在 rho_ref 未核实前不填数值 L、不输出结构。密度 ±3% 是有意环境覆盖，不表示三点均为真实平衡液体密度。

## 首选采样器：原始 MACE-MP-0b3-medium

来源与许可已取得官方证据：
- MACE-MP 官方 main 查询得到提交 `16a9f178706ce053f3ca8531efbab00a305d0854`。
- 固定提交模型目录： https://raw.githubusercontent.com/ACEsuit/mace-mp/16a9f178706ce053f3ca8531efbab00a305d0854/README.md ，HTTP 200，SHA256 `3981e18e5b2d83bd6223ffc0724a261c638a9a0c812eddad91645821c91a3a12`。表列 MACE-MP-0b3 覆盖 89 元素、MPTrj/PBE+U、要求 >=0.3.9、模型许可 MIT。
- 固定提交许可： https://raw.githubusercontent.com/ACEsuit/mace-mp/16a9f178706ce053f3ca8531efbab00a305d0854/LICENSE ，HTTP 200，SHA256 `31ea0ccf7bc19797081bff51c7eff3a3927c8cb6c1d7a726dec3d157b436da1c`。MIT 要求复制/分发软件或其主要部分时保留版权和许可声明。
- MACE v0.3.16 官方 loader： https://raw.githubusercontent.com/ACEsuit/mace/v0.3.16/mace/calculators/foundations_models.py ，HTTP 200，`medium-0b3` 指向 https://github.com/ACEsuit/mace-mp/releases/download/mace_mp_0b3/mace-mp-0b3-medium.model 。
- MACE v0.3.16 代码许可： https://raw.githubusercontent.com/ACEsuit/mace/v0.3.16/LICENSE.md ，HTTP 200，MIT，SHA256 `42137790f854ae2b9d29a0a72a4da3f6fb9a21d820b3f276f23b0575af72e86e`。
- 项目现有模型文件 `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model` 本次只做 SHA256，得到 `2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012`，与归档一致；未修改/推理/再次下载模型。远程 release 资产的独立散列本次未下载复算。

推荐：原始模型作为未来独立批准的**结构采样器首选候选**，复用已核验权重和 MACE 0.3.16 环境，不把 V18 的力退化模型用作默认液态生成器。官方训练集与 PBE+U 信息只说明来源，不证明银液态采样物理准确，不意味着形成有效 DFT 标签。

## 备选采样器：Foiles/Baskes/Daw Ag EAM

具体文件源可核对，已从 LAMMPS 官方仓库取得：
- LAMMPS develop 查询得到提交 `c71ef819ab3785bba81c31c0ef5fce84603608d3`。
- https://raw.githubusercontent.com/lammps/lammps/c71ef819ab3785bba81c31c0ef5fce84603608d3/potentials/Ag_u3.eam ，HTTP 200，36549 bytes，SHA256 `991e5f6a33c24908227d9c749b6ad9a8266200ce2955335ec113e9dd06116777`。
- 文件首行：`DATE: 2007-06-11 UNITS: metal CONTRIBUTOR: Stephen Foiles ... CITATION: Foiles et al, Phys Rev B, 33, 7983 (1986) COMMENT: Ag functions (universal 3)`；47 Ag，DYNAMO funcfl 单元素 EAM，库文件质量值约107.87。
- https://raw.githubusercontent.com/lammps/lammps/c71ef819ab3785bba81c31c0ef5fce84603608d3/LICENSE ，HTTP 200，GPL version 2，SHA256 `be38e38d9482c2beae35408e413045b4363e0e067b0d9f31549047948637a7ed`。
- https://raw.githubusercontent.com/lammps/lammps/c71ef819ab3785bba81c31c0ef5fce84603608d3/potentials/README ，HTTP 200，SHA256 `f672e30ad774c07c4c0a858d7013159472f4638f321225c8262e96a926fa413f`。官方提醒势文件主要展示类型、文件存在不保证转录和适用性，应与发表结果核验。

许可结论：已识别官方仓库 GPL-2.0 许可；Ag_u3.eam 本文件未包含单独的作者许可声明，不能声称它已取得独立 MIT/无条件公共域许可。如采用/分发这个具体候选，应保留官方仓库许可与来源，不绕过其条款；需要记录选用者对具体文件许可适用范围的复核。作为备选当前仍未完成液态适用性和数值实现核验。**本次未安装或运行 LAMMPS**；即使将来通过 ASE EAM 实现读取该文件，也必须独立授权采样、验证该实现的 EAM 格式/单位/插值与同一势的参考结果，不可静默换采样器。

## 液态采样判据对来源核验的影响

1250 K 接近实验熔点，64 原子体系从 fcc 直接加热至1250K可能保持亚稳固体。因此未来获批后应先通过批准的高温熔化/充分混合轨迹确认液态，再独立冷却/稳定至目标温度；高温目标、时间步长和过程预算需单独预登记。不能用静态 fcc 加随机位移并标记1250K代替液态。

两 seed 必须独立初始化完整轨迹，不能从同条轨迹复制终帧仅改seed；同初始化父族的改seed/密度/扰动仍需分组隔离，不自动创造独立测试来源。结构采样器的能量与力只作采样/稳定性诊断，全部不得作为DFT真值。参考rho、采样许可/来源、真实轨迹及液态QA仍须全部通过后才具有结构就绪条件。
