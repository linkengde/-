# 六个液态 Ag pilot：批准与放行门

pilot 六个结构是 `MASTER_TASK_MATRIX.csv` 的六行子集，不是另加六个结构任务。只有用户另行明确批准实际结构生成与六个 DFT 单点后，D/E/F 才能运行它们。当前全部坐标和 DFT 输入均 `STRUCTURE_NOT_READY`。

## 六个 pilot task ID

| 云端 | Task ID | 温度 | 相对密度 | seed | 结构状态 |
|---|---|---:|---:|---:|---|
| D | `C-D-AGL-1250K-rho097-seed01` | 1250 K | 0.97 × rho_liquid(1250 K) | 125009701 | `STRUCTURE_NOT_READY` |
| D | `C-D-AGL-1250K-rho103-seed02` | 1250 K | 1.03 × rho_liquid(1250 K) | 125010302 | `STRUCTURE_NOT_READY` |
| E | `C-E-AGL-1250K-rho097-seed02` | 1250 K | 0.97 × rho_liquid(1250 K) | 125009702 | `STRUCTURE_NOT_READY` |
| E | `C-E-AGL-1250K-rho100-seed01` | 1250 K | 1.00 × rho_liquid(1250 K) | 125010001 | `STRUCTURE_NOT_READY` |
| F | `C-F-AGL-1250K-rho100-seed02` | 1250 K | 1.00 × rho_liquid(1250 K) | 125010002 | `STRUCTURE_NOT_READY` |
| F | `C-F-AGL-1250K-rho103-seed01` | 1250 K | 1.03 × rho_liquid(1250 K) | 125010301 | `STRUCTURE_NOT_READY` |

每个为 64 Ag atoms、PBC=xyz、不同预定采样 seed/parent；实际轨迹、RNG 实现、来源哈希和坐标尚不存在。相对密度不是绝对 cell 体积。

## 放行前必须完成

1. **密度与来源：** 找到可复查的 Ag 液态 1250 K 密度表/方程，记录单位、插值和不确定度。已知约 1234.93 K 熔点值不能替代液体密度。旧有“熔点附近约 9.3 g/cm³”只是量级线索，当前不能拿它构造 DFT cell。
2. **sampler 与结构生成：** 指定生成器、权重/参数来源、版本、许可和 SHA。原始 MACE-MP-0b3-medium 或来源和许可经核验的 Ag-only EAM 可以提议候选坐标；其能量/力不能当 DFT 标签。结构生成 MD/AIMD 必须单独获得用户授权，本任务未批准它。
3. **真实液态证据：** 轨迹进入稳定取样区；温度、密度、能量、应力/压力稳定；RDF 与 S(k) 不保留长程 fcc Bragg 特征；保存 RDF 第一极小值处配位分布、Q6/CNA、非平台 MSD 和 seed 间 RDF/CN/SOAP 多样性。接受阈值在打开 DFT 结果前预注册。发生结晶的 1050 K future candidates 需标实情，不作为液体标签。
4. **身份、重复与角色：** 冻结 raw/candidate/trajectory/final extxyz SHA、frame index、元素和 atom ID 顺序、cell/PBC、seed/parent。精确哈希与周期近重复检查通过；全六个 parent 先登记为 proposed train-only 候选，不能分到 dev/test。现有 15 个 Ag 与 10 个 Ti₃SiC₂ 数值控制保持原角色。
5. **主机/方法：** D/E/F owner 只读报告活跃/排队/未启动任务、CPU/内存/盘、GPAW/ASE/PAW/MPI/ABI/许可。PBE/PAW1.2.1/PW500eV/σ0.10eV/SCF1e-5 固定；环境不兼容则停，不私自安装/换参数。独立运行区和保留空间就绪。
6. **k 网格：** 一个筛过液态条件的同一快照比较 Γ 与 2×2×2；如果 ΔE>2 meV/atom、force RMS>0.01 eV/Å 或 max force delta>0.03 eV/Å，则增至 3×3×3 或更密并复核。不能用不同 seed 结构代替固定结构 k 收敛对照。mesh 对照不产生新 parent/构型标签，须单列计算耗时与哈希。
7. **批准/预算：** C 收齐上述材料后提交具体 pilot 输入清单、source hash、总预算和停止条件；用户另行明确批准后才可开始。单例暂定上限 12 h 或 160 SCF iterations；首例超过 8 h 即停后续并复估，不自动重启。

## 每个结果的验收

所有 6 个标签都必须 SCF 收敛；native `energy`、force-consistent `free_energy` 和全部三维原子力有限；单位 eV / eV·Å⁻¹；输入坐标顺序、cell/PBC、参数、k 网格与冻结清单匹配；mesh 对照误差满足预算；source/input/result SHA-256 及归档 verifier 全 PASS。任何 NaN/Inf、SCF 未收敛、来源/license 不明、parent 重复、SHA 不符或预算超限，该标签隔离并停止本批剩余新任务。

C 审查完成后，pilot 通过也不自动放行剩余 24 项；用户必须进一步批准正式批次。
