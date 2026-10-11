# D/E/F 最新验收与 Ag64 首项暖机放行总结

更新：2026-10-11 UTC

## 节点汇总

| 节点 | owner 分支当前核验头 | 势 SHA / ASE 元数据解析 | 合成归档往返 | 当前预检 | C 放行 |
|---|---|---|---|---|---|
| D | `2c8b6f41ec8a7605ddc06f7218014c55957166ac` | PASS / PASS | PASS，20480-byte fixture；本轮报告已 fresh fetch 校验 | 旧 HOLD 配置下仅 C release 子项失败，其余通过 | 仅首 task 1000 步，需新 SHA 下本机 preflight PASS |
| E | `39734568216e5a3d9853d7d052d28f9799172598` | PASS / PASS | PASS，20480-byte 确定性 fixture，从全新仓库恢复并逐字节及 SHA 核对 | 旧 HOLD 配置及旧配置 recovery 状态下 HOLD | 不放行 |
| F | `2ddb39a899decd6b3cc8d6c002d48291c9a5e632` | PASS / PASS | PASS，既有 16415-byte fixture及本报告 fresh-fetch 核验 | 旧 HOLD 配置下仅 C release 子项失败，其余通过 | 不放行 |

三个 owner 回执的文件清单由 C 独立从远端提取并校验：D 报告/原始 preflight SHA 分别为 `6ce619eea8da27e3aff97a93a3893d09187cc50eb2a41410466a051072f8d2f7` / `dffd8945cdc7d61cb7c76811321965b4fe1b12ad6e7307e3c9b3e21305eb2cc7`；E 当前验收清单 7 项全通过；F 当前验收清单 6 项全通过。分支和逐文件证据见本目录 `OWNER_EVIDENCE_INDEX.csv`。

三台机器报告的 CPU quota 均为 4 核，cgroup 内存上限均为 32 GiB，工作盘余量约 28.79–28.87 GB，inode 余量约 203 万。工作区是临时 overlay；配置继续 `persistence_verified=false`，远端归档是必须的恢复副本。容器报告不能证明物理宿主机独占。

势函数 cutoff `5.995011 Å` 大于候选 64 原子立方盒的半边长。C 已在既有 `SAMPLER_SOURCE_REVIEW.md` 审核本执行器所用 ASE `PrimitiveNeighborList` 的整数周期 offset 和多周期像行为：该邻居表支持 cutoff 跨半盒长，且执行器保留同一原子的不同周期像。因此这里不构成采样实现阻塞；不可替换成 minimum-image 或按 atom index 去重。这个实现审查不代表有限尺寸熔点、扩散和凝固已通过物性验证。

## 当前配置

- 现行配置：`CODEX_WARMUP_RELEASE_CONFIG.json`
- SHA-256：`4aee9d64859736f35459259625598f767d5f3e0d8b33bc71f0abb5baf31639e2`
- 只授权 D 的 `C-D-AGL-1250K-rho097-seed01`，暖机最多 1000 步、0.5 fs。
- D 的 C host release 为 true；E/F 为 false；D 第二任务、E/F 四任务均 HOLD。
- continuation 为 false：不自动启动余下 119000 步，也不授权 DFT、MACE 训练、TTM 或 LAMMPS。

## D 唯一下一步

D 必须先按 `C_RELEASE_ATTESTATION_D_WARMUP_20261011.md` 中的命令，以此完整 config SHA 和精确本机势路径重跑 `node_preflight_codex.py`。仅当当前会话输出 PASS 时，D 才可运行冻结的 `run_sampling_codex.py`。运行器会再次核验配置/代码 SHA、任务身份、当前 owner 分支 head、本机势、工作根/归档根及资源条件，并在完成后自动执行归档与 fresh-fetch 恢复校验。分支 HEAD 已变化、任一校验非 PASS、预算超限或备份失败均须停止并报告，不调整参数绕过。

本次 C 仅签署配置，没有远程访问或启动 D 上的过程；当前没有新的 Ag 原子坐标、轨迹、checkpoint 或模拟结果。
