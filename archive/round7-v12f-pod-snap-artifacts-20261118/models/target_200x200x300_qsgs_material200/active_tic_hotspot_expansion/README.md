# Ti–C 骨架热点标签候选

本目录记录全模型冷态稳定性筛查发现的 Ti–C 短接触，保留 Si01 为开发留出，且不改写 v12f/v12g。

## 发现

- 配对原子 ID 62864 (Ti) / 62871 (C) 在完整初始结构中的间距为 2.0793 Å。
- v12g 完整模型无热源 10 fs 后，该对间距降至 1.6826 Å。
- v12f 同一筛查有 47 对 Ti–C 小于 1.8 Å；v12g 有 53 对。最低邻距分别约 1.6865 Å 和 1.6826 Å。
- 对上述初始 Ti–C 对，v12g 的全模型净力给出沿 Ti→C 方向的质量加权相对加速度投影为负，意味着该瞬时力倾向于缩短这对原子间距。此值是基于原子合力的诊断量，不是分解后的成键力；需要 DFT 标签核实。

## 构型与标签

`candidates/tic_pair_initial.extxyz` 和 `candidates/tic_pair_cold_10fs.extxyz` 都是 36 原子、5.5 Å 截取半径、22 Å 真空立方盒的局部构型；两中心原子都标为 DFT 力目标。先计算初始态，再计算 10 fs 压缩态。采用 GPAW PBE/PW350/Gamma、0.10 eV Fermi–Dirac 展宽，与 Si 热点标签保持一致。标签收敛后再形成独立 v12h 训练候选，局部非目标原子沿用 POD v10 力锚点，并重新做 FitPOD 列顺序审计。

`candidate_manifest.json` 和 `../cold_v12f_initial_458k_1fs/short_contact_pair_audit.json`、`../cold_v12g_initial_458k_1fs/short_contact_pair_audit.json` 记录来源、映射和邻距。Si01 ID 104359 不进入训练集；原 00/02/03 训练角色不变。
