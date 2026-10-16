# H 修订：48 原子 Ti₃SiC₂ 静态小位移 DFT 单点预检

**日期：** 2026-10-16

**范围：** 单个已审核候选的离线预检与资源审计；未启动 DFT、MD、AIMD、TTM、松弛或训练。

**最终状态：** `DATA_ROLE=HOLD`；`RESOURCE_READY=NO`；`DFT_EXECUTION_READY=NO`；`execution_authorized=false`。

## 审核来源和候选绑定

已从 GitHub 读取并核验 C 提交 `25bb74065797d516f1e8c8bbf9b68a4b558aefef`，其审核报告 SHA-256 为 `1c6a4ff44882b5553048dd009a5fc6c0ff5ec8939210fdfae12f0004fcc7ed72`，决策 JSON SHA-256 为 `f7e71dd54c6590f5ff2612099489b0c54ce44d1384ee13d61f45762c86c0e832`，审核包 SHA 清单文件 SHA-256 为 `aa1039df5642b28560aa1b8ae25aec0ae13045d5c148cf3620762184b18692e7`。C 结论为 `C_SCIENTIFIC_REVIEW=NEEDS_REVISION`、`DATA_ROLE=HOLD`、`KGRID_CANDIDATE=REVISE`、`RESOURCE_READY=NO`、`DFT_EXECUTION_READY=NO`，没有签发执行授权。

C 审核绑定的 H 基准提交为 `e923008453d71077421b651cf19c6705420d2356`。单一候选文件为 `research/mace-v12-transfer/coordination/reports/chatgpt-nep-ace-dft-repair-20261016/NEW_DFT_GEOMETRY_CANDIDATES/structures/TSC48_static_perturb_std0p025_seed2026101603_UNLABELED.extxyz`，生成源提交 `cc65a3e6b75be28c06fda8608b4064662d714117`，SHA-256 `da3820c477023e0bb423f87550ed82ae0c2a571f53d294c055f78d7d052d5e8b`，大小 2,572 字节。它有 48 个原子（Ti24Si8C16），晶胞和 PBC 与母体一致，PBC=T T T；母体 SHA 为 `157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4`。C 已核实几何基本检查通过，最近周期 Ti–C 距离 2.0016108 Å。候选没有 DFT 标签，属于 `TSC_BULK_PARENT_OR_CONTROL_FAMILY`；保持 HOLD，不分配 TRAIN/DEV/TEST。

原 H 报告 SHA-256 为 `eec7bb2b0af188a9cd34f8199c234ccca361a53af1e533a23d7f06c501d31aa8`。本次报告作为该报告的增补修订，不覆盖原文、输入结构或历史 DFT 输出。

## k 点定义修正

| 项目 | 原 H 报告 | 修订结果 |
|---|---|---|
| 网格 | 将偶数 10×10×4 错标为 Γ 居中 | 保留尺寸，采用母体实际使用的标准 Monkhorst–Pack 网格 |
| GPAW API | — | 母体 runner 使用 `kpts=tuple(record['method']['kpts'])`；没有 `gamma=True` 参数，也没有额外 offset |
| GPAW 日志 | 原解释不正确 | `MonkhorstPack(size=[10,10,4])`、shift `[0,0,0]`；偶数网格不含 Γ 点 |
| 坐标证据 | 未给坐标表 | 离线生成完整 400 点坐标；规范化坐标列表 SHA-256 `530f266922e6154492fa0cf6608bf4c1b87e5e2b3a03f7176287c5c231bcaaec`；母体日志展示的 11 行均与该 full-BZ 网格一致 |

详细坐标与解析校验见 [KPOINT_MESH_AUDIT_10x10x4.json](KPOINT_MESH_AUDIT_10x10x4.json)。母体日志只显示 200 个 IBZ 点中的 11 行，因此完整 IBZ 行逐项比对仍受日志省略限制；网格类型、大小、零 shift、runner 参数及所有可见坐标均相互吻合。

同几何母体 8×8×4→10×10×4 的能量与力差均通过原定阈值，支持把 10×10×4 作为该候选的初始网格密度。**母体网格收敛通过不等于当前扰动几何已完成 k 点收敛验证**；C 也将 candidate k-point convergence 标记为 `UNVERIFIED`。

## 计算配置锁

候选配置 SHA-256：`e3edf53ebc1fcf1b3bf7383a762d6fa29e509cb2cfa307cb9d0bf136e10306fb`。沿用母体已归档方法：PBE/PAW、PW 500 eV、Fermi–Dirac 0.10 eV、非自旋极化、energy/density/eigenstates 阈值各 `1e-5`、PPCG `niter=2`、Pulay `beta=0.05, nmaxold=8, weight=100`、最大 SCF 迭代 160；原子坐标和晶胞固定。GPAW 26.7.0、ASE 3.29.0、gpaw-data 1.2.1、OpenMPI 5.0.7。Ti/Si/C PAW SHA：`bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f` / `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0` / `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983`。当前 runtime 脚本 SHA：`e32d19eef429939bcb3d0aa2d365d063b6dca16b4b4052b02da226fd7247e43b`。

执行资源建议仍为 MPI=4、OMP=1，但本 Codex 会话的 4-rank MPI 初始化在 PMIx listener socket 处失败（errno=1，exit 213）。C 归档的历史四 rank 测试确实通过，rank 0–3 均报告 size=4；该历史证据不代表当前 sandbox 的启动测试成功，也不显示 ABI 不兼容。当前 MPI 一致性状态为 `UNVERIFIED/BLOCKED`。

## 当前资源与恢复准备

本次当前环境快照见 [RESOURCE_PREFLIGHT_REVISION.json](RESOURCE_PREFLIGHT_REVISION.json)：CPU quota 4 核、cgroup 内存上限 32 GiB；`/workspace` 可用 18,775,195,648 字节。工作区为任务隔离 overlay，`/mnt/e-work` 不存在，`DFT_WORK_ROOT` 未设置，因此没有会话结束后的持久化依据。C 的旧快照曾报告可用磁盘 28,828,061,696 字节；本次现场可用空间更少。30 GB 只是建议余量，不是 GPAW 的硬性最低值，但当前不足以证明本单点的 scratch/checkpoint 写入预算安全。

用母体日志中的 200 个 IBZ 点、235 bands、14,525 个平面波系数和 complex128 计算，完整波函数负载的稠密数据量级代理为 10,922,800,000 字节（10.17 GiB）。这不是实测 GPW 文件大小，也不包括元数据、scratch 和写入临时占用。单份代理量小于当前空闲字节，但不证明实际 checkpoint 及写入余量足够。母体参考 runner 每 20 步将 `mode='all'` 写到同一个本地 `state.gpw`；它没有经验证的 restart 分支，也未证明写失败时旧 checkpoint 一定完整。当前没有持久卷、可恢复 checkpoint、外部墙钟上限或实测 RSS 证据。

**`RESOURCE_READY=NO`**：需由 C 审查单候选磁盘例外或持久存储方案，并冻结 scratch/检查点容量、checkpoint 持久性、验证恢复方式与外部 walltime 证据。不能仅用旧母体的 3.26 小时运行时长代替作业墙钟配额。

## 复审请求和停止状态

请 C 复审：标准 MP 定义及 k 点坐标证据；当前环境与历史四 rank 证据的差异；单候选磁盘预算、持久性和恢复路径。候选属于母体父族；TRAIN/DEV/TEST 继续 HOLD。

在 C 复审、资源放行及用户另行明确授权前，禁止启动该单点。此报告只完成修订和预检，`DFT_EXECUTION_READY=NO`。
