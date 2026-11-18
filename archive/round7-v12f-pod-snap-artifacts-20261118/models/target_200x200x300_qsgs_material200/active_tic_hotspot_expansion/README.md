# Ti–C 骨架热点标签候选

本目录记录全模型冷态稳定性筛查发现的 Ti–C 短接触，保留 Si01 为开发留出，且不改写 v12f/v12g。

## 发现

- 配对原子 ID 62864 (Ti) / 62871 (C) 在完整初始结构中的间距为 2.0793 Å。
- v12g 完整模型无热源 10 fs 后，该对间距降至 1.6826 Å。
- v12f 同一筛查有 47 对 Ti–C 小于 1.8 Å；v12g 有 53 对。最低邻距分别约 1.6865 Å 和 1.6826 Å。
- 对上述初始 Ti–C 对，v12g 的全模型净力给出沿 Ti→C 方向的质量加权相对加速度投影为负，意味着该瞬时力倾向于缩短这对原子间距。此值是基于原子合力的诊断量，不是分解后的成键力；需要 DFT 标签核实。
- 另抽取独立位置的 Ti ID 324028 / C ID 324034；v12g 冷态 10 fs 后间距为 1.6868 Å，局部截取 38 个原子，用于检查失稳是否跨区域重复。

## 构型与标签

`candidates/tic_pair_initial.extxyz` 和 `candidates/tic_pair_cold_10fs.extxyz` 是 36 原子局部构型；`candidates/tic_pair2_cold_10fs.extxyz` 是第二处 38 原子局部构型。均使用 5.5 Å 截取半径和 22 Å 真空立方盒，两中心原子标为 DFT 力目标。DFT 使用 GPAW PBE/PW350/Gamma、0.10 eV Fermi–Dirac 展宽，与 Si 热点标签保持一致。初始态和主热点压缩态均已收敛；独立位置的第三个标签仍在计算。主热点压缩态 1.683 Å 时，DFT 的质量加权相对净力投影为 +0.625，v12g 为 -0.0106，说明 v12g 对这个压缩接触给出的分离趋势明显不足。该投影基于原子合力，不是成键力分解。所有标签完成后再拟合隔离的 v12h 候选；局部非目标原子沿用 POD v10 力锚点，并审计 FitPOD 力列顺序。

`candidate_manifest.json` 和 `../cold_v12f_initial_458k_1fs/short_contact_pair_audit.json`、`../cold_v12g_initial_458k_1fs/short_contact_pair_audit.json` 记录来源、映射和邻距。Si01 ID 104359 不进入训练集；原 00/02/03 训练角色不变。
