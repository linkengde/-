# D 首个 Ag64 暖机：C 归档与执行审查

审查日期：2026-10-11 UTC  
范围：仅读审 D owner 分支的暖机归档、readback 收据与记录；未启动或续算 MD、DFT、训练或其他计算。

## 决定

**暖机操作完成：PASS。远端归档与恢复：PASS。结构/液态科学 QA：未通过（未执行），液态身份仍未确认。** 这次 1000 步任务满足其有限的暖机和归档检查，可作为运行流程/资源试跑记录。它不能当作平衡液体、正式 Ag 构型或 DFT 候选；当前不放行续算或其他五条暖机。

温度日志存在需要在后续方案中处理的启动瞬态：热段 target 为 1600 K，step 80（40 fs）记录峰值 4266.30 K；step 1000（500 fs）仍为 1855.25 K，约高于 target 15.95%。签署配置没有温度容差/平衡通过阈值，因此此结果不构成执行器违约；但它不支持“温度已稳定”的说法。应先审定温度/热化 QA 门槛，再决定其余暖机是否继续。

## 远端归档和恢复

D 分支 `cloud-d-b1-dft` 最新头为 `145ee8050cf95a3a1988e868c6b935e10538f723`。它是只读核验的目标，不由 C 改写。暖机归档先在 commit `11803909226538ab83847d1f06827a8292501139` 写入；readback 回执在 `145ee8050cf95a3a1988e868c6b935e10538f723`。

C 另用全新 clone/fetch 取回该分支，再独立检验归档：

- 归档 `warmup-b62a818e58a2cf42.tar.gz`：100,700 bytes，SHA-256 `b62a818e58a2cf423e4cb1b9ed0468ae73ad578a75c4b2cb9fb3f28baf75e7b1`。
- `archive-manifest.json`：SHA-256 `ab9f0dffa8ef145eac7792dfe1e14e07c1354df1c2152dafec97bad1dc3a57b1`。
- `REMOTE_READBACK_RECEIPT.json`：SHA-256 `8bc2b17738f5e214098964c0635ef3505134726f8b4d8e7323cd4c5b03591bd4`。
- 外层归档 SHA、manifest SHA 和 18 个成员的路径、大小及逐成员 SHA 全部复核通过；readback 收据记录 fresh-fetch restore 18/18。
- checkpoint `state_arrays.npz`、单帧 `state_atoms.traj` 与总轨迹末帧的 positions、cell、numbers、masses、PBC、atom IDs 和 momenta 一致；坐标/速度数组有限值检查通过。
- 归档不包含 `Ag.eam.alloy` 或论文文件；manifest 标记 `potential_file_included=false`、`unlicensed_documents_included=false`。

格式小项：归档提供完整的 JSON `archive-manifest.json`，不是名为 `SHA256SUMS.txt` 的清单；tar 内也没有独立 `README.md`。这不影响本次 18 成员恢复和 SHA 验证，但与统一归档协议要求的显式 README 文件不完全一致，应在后续归档实现中补齐。

## 任务身份、参数和运行资源

归档身份与本次授权一致：D / `C-D-AGL-1250K-rho097-seed01`，配置 SHA `2e122b258a48fe54bf603709b5505a301c1257f398c5625eecc2ca94998c0052`，Ag 势 SHA `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`，64 个 Ag、10.826824658074 Å 立方盒、三轴 PBC、密度 9.03273341876 g/cm³、seed 125009701、步长 0.5 fs。轨迹时间为 0 至 500 fs，共 11 帧（每 100 步）；热日志 51 条（每 20 步）。运行完成 1000 步，0.5 ps。归档记录中的任务矩阵 SHA 与原六项授权清单一致，归档记录的四个科学/归档执行代码 SHA 与 C 签署配置一致。

暖机 wall time 200.307 s，约 4.992 steps/s；记录 RSS 102,563,840 bytes，结束时空闲盘 28,792,524,800 bytes。它们都在配置的 3,600 s、8 GiB RSS 和 2 GiB 空闲盘线内。归档与恢复也低于 50 MiB 限制。运行记录为 `DFT_STARTED=false`，source manifest 为 `DFT_LABELS=ABSENT`。

## 静态结构与热日志检查

C 对归档轨迹做了只读数据检查：11 帧均为 64 Ag、完整 PBC，原子坐标和晶胞均有限；全轨迹最小周期像原子间距为 2.200686 Å，末帧为 2.335474 Å，未低于执行器 1.8 Å runtime minimum。它只说明这些已记录帧没有出现该短程距离门禁违例，不替代完整结构 QA。

`sampler_state.json` 明确保留 `QA_status=NOT_PERFORMED_BY_SAMPLER`、`STRUCTURE_PASS=false`；`source_manifest.json` 也标记 `liquid_QA=NOT_PERFORMED` 和 `DFT_LABELS=ABSENT`。本审查没有以 1250 K 或单段轨迹认定液态，也没有计算 RDF/Q6/MSD 来把暖机升级为液态验收。

## 下一步边界

- D 首项 1000 步暖机及其归档/恢复审查结案为操作流程 PASS。
- 该构型仍是 `liquid_pilot_candidate_not_training_or_independent_test`；未获 DFT 准入。
- D 的剩余 119,000 步、D 第二项、E/F 四项和其他节点任务均未由本报告放行。
- 在决定后续任务前，需先设定并审查高温瞬态/热化 QA 规则，并为后续科学归档补入独立 README 文件。
- A/B 工作区和数据未触碰；G/H 仍待命。
