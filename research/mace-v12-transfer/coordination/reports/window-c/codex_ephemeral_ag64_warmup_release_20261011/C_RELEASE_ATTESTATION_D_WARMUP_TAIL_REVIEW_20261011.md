# C 补充签署：Ag setfl 径向尾部审查与 D 首项暖机配置修订

日期：2026-10-11 UTC
决定：`TAIL_POLICY_APPROVED_FOR_EXACT_FILE`; 仅保留此前 D 首项暖机范围。

## 审查结论

D 的执行器在第 0 步因 `PUBLISHED_SETFL_RADIAL_TAIL_C_REVIEW_REQUIRED` 安全停止。根据 D 提供的执行回执，没有生成原子坐标、采样目录、checkpoint 或归档，也没有启动 MD/DFT。门禁没有被绕过。

C 已复核仓库中既有的 `SAMPLER_SOURCE_REVIEW.md`、`SAMPLER_PARAMETER_VERIFICATION.json` 和 `SAMPLER_INSPECTION.json`。这些记录的 SHA-256 分别为：

- `SAMPLER_SOURCE_REVIEW.md`: `8883a516a9227fa7c4681135fd707ca97c525285285ed1a3a59b25370c7eacfa`
- `SAMPLER_PARAMETER_VERIFICATION.json`: `738c2c5ee1bb0df07d72e308c5b5caab0c1f81c00b0ca30843208277e3c9d930`
- `SAMPLER_INSPECTION.json`: `0f9cc608d026cb36a801d282181f39473f761f4d3bcf42b1f1c2ee0b48c6db1c`

它们对应的唯一势文件是用户提供的 `Ag.eam.alloy`：780452 bytes，SHA-256 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`，setfl/alloy，ASE 3.29.0 可解析。C 在当前副本重算 SHA 一致，并以 ASE 3.29.0 仅解析表格和样条元数据；没有创建 Atoms、计算能量/力或运行 MD。

文件给出的 `Nr=10000`、`dr=0.0005995011000293093 Å` 和 cutoff `5.995011000293092 Å` 满足 cutoff=`Nr×dr`。ASE 解析的最后径向表点为 `(Nr−1)×dr=5.994411499193063 Å`，因此表格终点到 header cutoff 的间距恰为一个 `dr`（仅有浮点舍入差）。沿 ASE 3.29.0 对该原样表格的 cubic-spline 插值，此窄区间 11 点检查得到 `max|electron density|=1.979072636197065e−14`、`max|pair phi|=2.0656530567331374e−13 eV`；端点值分别约 `−1.40e−16` 与 `1.46e−15 eV`。尾段导数有限且接近数值零，但不声称函数或导数数学上严格为零。

因此，已发布文件的 header cutoff 比最后径向网格点多一个网格间隔，且本执行器采用的 ASE 插值在这个确切文件的间隔上数值接近零。ASE 3.29.0 的实现使用 SciPy `InterpolatedUnivariateSpline(k=3)`。11 点内插检查中最大电子密度导数绝对值为 `6.764210394245453e−11`、pair phi 导数绝对值为 `7.059778089168755e−10 eV/Å`。C 批准窄策略 `C_APPROVED_PUBLISHED_SETFL_ONE_BIN_ZERO_TAIL`，**仅适用于上述确切 SHA、`alloy` 格式、ASE 3.29.0、且通过运行器已有的 gap/尾段阈值检查**。这不是对其他 setfl/EAM 文件、任意外推区间或其他程序插值行为的通用认可。没有改写势文件、PLT、cutoff、势函数值、采样器或科学参数；文件的再分发状态仍为 `UNVERIFIED_DO_NOT_REDISTRIBUTE_POTENTIAL_FILE`。

运行器已经实现显式策略校验：必须是 `alloy` 格式、尾段不超过一个 `dr`、表函数尾值低于原有 `1e−10` 阈值，否则仍然 HOLD。D 的停止原因是签署配置缺少已审策略字段，而不是尾值超限。现行配置现补入该策略及证据 SHA；没有放宽运行器判据。

## 修订配置和范围

- 新配置：`CODEX_WARMUP_RELEASE_CONFIG.json`
- 新 SHA-256：`2e122b258a48fe54bf603709b5505a301c1257f398c5625eecc2ca94998c0052`
- 上一份已签配置原样保存在 `CODEX_WARMUP_RELEASE_CONFIG_PRE_TAIL_POLICY_20261011.json`，SHA-256 `4aee9d64859736f35459259625598f767d5f3e0d8b33bc71f0abb5baf31639e2`

本签署由 C 审查记录、配置 SHA 和本分支不可变 Git 提交固定；它不是公钥数字签名。

新配置只补充上述势文件尾部策略、证据引用及签署记录。仍然只授权 `C-D-AGL-1250K-rho097-seed01` 的 `warmup` 阶段，最多 1000 步、冻结步长 0.5 fs。`persistence_verified=false` 保持不变；D 第二项、E/F、其余五个任务、完整 119000 步续算、DFT、训练、LAMMPS 与 TTM 仍为 HOLD。

D 上一次 preflight PASS 绑定旧配置，不能替代新配置下的当前资源检查。D 必须先从 C 分支取得新配置及同版执行包，在本节点同一会话重新运行只读预检：

```bash
python node_preflight_codex.py \
  --config CODEX_WARMUP_RELEASE_CONFIG.json \
  --config-sha256 2e122b258a48fe54bf603709b5505a301c1257f398c5625eecc2ca94998c0052 \
  --owner D \
  --task-id C-D-AGL-1250K-rho097-seed01 \
  --potential-path /workspace/WMH2006_Ag-local-validation/Ag.eam.alloy
```

只有当前结果为 `PREFLIGHT_STATUS=PASS` 且门槛没有变化时，才可在 D 自己的 Codex 会话运行原执行命令：

```bash
python run_sampling_codex.py \
  --config CODEX_WARMUP_RELEASE_CONFIG.json \
  --config-sha256 2e122b258a48fe54bf603709b5505a301c1257f398c5625eecc2ca94998c0052 \
  --owner D \
  --task-id C-D-AGL-1250K-rho097-seed01 \
  --stage warmup
```

运行器仍会二次检查势 SHA、策略、代码 SHA、远端 D 分支头、资源及归档门。任一检查不是 PASS、远端 head 已改变、输出超预算或备份/恢复失败，都必须停止并报告，不得通过改势、改 cutoff 或弱化检查继续。本补充签署没有远程启动 D，也没有产生新的采样结果。
