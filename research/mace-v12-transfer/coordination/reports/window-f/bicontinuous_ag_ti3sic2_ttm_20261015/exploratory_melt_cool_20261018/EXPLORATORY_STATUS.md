# Ag–Ti₃SiC₂方案A：短松弛与300 K探索性筛查

## 结论

模型副本已完成一次固定晶胞的有界 FIRE 短松弛，并启动了 300 K NVT 基线筛查。短松弛后的 Ti₃SiC₂ 几何骨架保持连通，主体 Ag 网络仍占绝大多数；但 300 K 前 0.2 ps 出现重复的局部极端受力，最大力约 66–67 eV/Å。故本轮在探索性最低安全门槛处停止，**未启动局部加热、冷却、TTM、DFT 或训练**。目前不能把这一模型和 calibration-v0 势函数组合用于熔池/凝固结论。

## 模型、势函数与隔离

- 输入为39,486原子方案A：Ti 13,576、Si 4,768、C 8,800、Ag 12,342；晶胞110 × 98.085134 × 100.920749 Å，`boundary f p p`。
- 结构副本与上传源模型 SHA-256 一致：`e2b62a68c46e23037370a6ca20951aa0774edea2d728d3c3aadbcba4b3755674`。原始模型文件未改写。
- 保持 calibration-v0 映射原样：Ag–Ag 使用 Ag EAM；Ti–Si–C 使用 TSCBOP Tersoff；Ag–Ti/Si/C 使用三组 Morse 交叉项。未与 MACE、EAM 或 TSCBOP 额外叠加。
- 势函数文件哈希与归档记录一致；calibration-v0 仍是未验证的探索性跨相势，不因短松弛能收敛或平均力降低而视为已验证。

## 短松弛

使用源模型独立副本、固定晶胞、FIRE、`dmax=0.05 Å`、最多500步/5000次力评估，目标 `Fmax≤0.5 eV/Å`。计算到步数上限而停止，未达到目标：最终 `Fmax=2.4501`、`P95=0.1619`、`P99=0.3836 eV/Å`，4个原子高于2 eV/Å，无原子高于5 eV/Å。

位移主要发生在Ag：Ag位移中位数1.297 Å、P99为4.218 Å，1,444个Ag原子移动超过3 Å。Ti/Si/C位移中位数约0.35–0.37 Å，P95均小于0.96 Å；只有4个Si移动超过3 Å。TSC图在3.30 Å连边阈值下仍为一个分量；Ag在3.05 Å下最大分量占99.708%。

同时，Ti–C最短距离由2.092 Å缩至1.729 Å；松弛后有20对Ti–C短于1.8 Å、1,668对短于2.0 Å。最短的两对（2047–2267、4409–4629）属于原始界面候选原子。Ag–C最短距离为1.906 Å（324对短于2.0 Å），Ag–Ti最短距离为2.027 Å。几何图连通不排除这些局部压缩。

## 300 K基线筛查

从 `post_relax.data` 副本开始，NVT目标300 K、时间步长1 fs。达到0.2 ps时中断，未完成平衡，也没有加热源。0.1 ps和0.2 ps的整体温度分别为199 K和263 K；Fmax分别为66.491和67.010 eV/Å，F>5 eV/Å的原子数分别为5和20。没有NaN或丢失原子的记录，但强受力反复出现，故不能把该分支计为通过的300 K稳定性测试。

典型高力点为C ID 19219：Fmax在最小化步475曾为约50.83 eV/Å，300 K筛查0.1 ps又达到66.49 eV/Å。该原子的初始最近Ag距离2.851 Å，属于相间界面候选区；0.1 ps时最近邻包含Ti ID 25627，距离1.809 Å，Ti ID 25628距离1.914 Å，以及Ag ID 38116距离1.981 Å。0.2 ps时局部C/Ti位点仍有约67 eV/Å的力。该重复局部响应是停止热驱动的直接理由，并非全模型骨架已经崩塌。

TSC在0.1/0.2 ps的2.70 Å几何图分别出现3/2个分量，最大分量分别含27,142/27,143个（总数27,144）；Ag在3.05 Å下最大分量仍分别占99.692%/99.643%。小分量随近邻阈值和热振动改变，不能据此声称完整化学键连通。

## 势函数与科学结论边界

E窗口统一静态对照提交 `dbcf2f7d6bf966b8fc0476224d2cba52c3866de0` 及本仓库已归档的 MACE 对照显示不同模型间存在残余力分歧；那些分歧不是当前结构的 DFT 标签。39,486原子结构没有DFT力标签，因此这里不报告DFT RMSE，也不以更低残余力推断模型更准确。当前异常可能涉及初始界面构型、TSCBOP相内响应及Ag–TSC Morse交叉项，离线结果本身不能唯一归因。

## 状态标签

- `SOURCE_MODEL_UNCHANGED=YES`
- `SHORT_RELAXATION_COMPLETED=YES; FMAX_TARGET_NOT_REACHED`
- `TSC_BACKBONE_GEOMETRY_CONNECTED=YES at 3.30 Å after minimization; 300 K graph threshold is vibration-sensitive`
- `AG_NETWORK_LARGEST_COMPONENT=99.708% after minimization at 3.05 Å`
- `300K_BASELINE=STOPPED_AFTER_0.2_PS; NOT_PASSED`
- `INTERFACE_FORCE_RISK=HIGH_UNRESOLVED`
- `READY_FOR_LOCAL_HEATING_COOLING=NO`
- `READY_FOR_STABLE_MELT_POOL_OR_SOLIDIFICATION_CLAIMS=NO`

## 文件

- `MINIMIZATION_OVITO.extxyz`：FIRE松弛的多帧轨迹，OVITO可读。
- `TRAJECTORY_300K_OVITO.extxyz`、`trajectory_300K.dump`：300 K部分轨迹；只有0、100、200 fs三个存档帧，不是完成的稳定性结果。
- `post_relax.data`、`post_relax.dump`、`post_relax.restart`：独立松弛副本输出；未覆盖源结构。
- `SCREEN_300K_AUDIT.md`、`RELAXATION_LOCAL_GEOMETRY_AUDIT.md`：摘要审计。
- `MINIMIZATION_SHORT_CONTACTS.csv`、`FORCE_SPIKE_LOCAL_NEIGHBORS.csv`、`300K_HIGH_FORCE_ATOMS.csv`：逐原子/近邻明细。
- `MINIMIZATION_STEP_SUMMARY.csv`、`300K_FRAME_AUDIT.csv`、`CONNECTIVITY_CUTOFF_SENSITIVITY.csv`：统计与阈值敏感性。
- `MODEL_SHA256SUMS.txt`：本分支文件清单与哈希。
