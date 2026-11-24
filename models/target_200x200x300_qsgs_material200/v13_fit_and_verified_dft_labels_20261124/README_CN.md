# V13 完整 DFT 标签与 POD 首轮拟合归档

归档日期：2026-11-24。包含 27 个通过审计的训练 DFT 标签、6 个独立来源验证标签和首轮四元素 POD 系数。

## 内容

- `training/`：27 个训练候选的几何、GPAW SCF 日志、标签 extxyz、标签 JSON。
- `validation/`：6 个独立父来源验证候选，严格保留为验证用途。
- `*_candidate_manifest_v13.json`、`dft_label_status_v13.json`：来源追踪与标签验收状态。
- `fit_v13_weight100/`：当前保留的 POD v13 描述符、系数、精确拟合输入、日志和逐构型误差。
- `fit_v13_weight300_comparison/`：能量权重 300 的模型选择对照结果。
- `V13_FIT_REPORT_CN.md`：中文结果与适用范围。
- `SHA256SUMS.txt`：本目录文件校验和（不包含该校验文件自身）。

## 审计与拟合摘要

标签验收：训练 27/27，验证 6/6。两组按父来源及规范几何哈希隔离。标签协议为 GPAW/PBE、PW350 eV、电子展宽 0.10 eV；该数值设置属于探索性统一协议，未宣称全面数值收敛。

权重 100 候选：训练能量/力 RMSE 0.068885 eV/atom、0.614599 eV/Å；验证能量/力 RMSE 0.105894 eV/atom、0.603560 eV/Å。能量验证误差略高于约 0.1 eV/atom 目标，Ti₃SiC₂ 热振动验证帧误差为 0.227383 eV/atom。验证池用于模型选择，不能称为最终盲测。权重 300 对照未改善验证结果。

本归档中的 V13 拟合候选未通过小体系无热源 NVE：液态 Ag、Ti₃SiC₂、Ag/Ti₃SiC₂ 界面和 Ag-rich Si 均出现塌陷或丢原子。当前不能用于长时加热或冷却。详细日志见 `cold_stability_screen/` 和 `V13_COLD_STABILITY_SCREEN_CN.md`。不要将“DFT 标签通过 SCF/哈希审计”解释为 DFT 数值设置普适收敛。
