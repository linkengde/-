# D 首项 Ag64 step 1000→2000 离线执行准备审查

## 决定

本轮完成只读源档核验、独立 step 边界门禁候选、合成故障测试、热化判据和预算。**未启动 MD/真实积分、DFT、训练或势文件推理。** `NEW_MD_RELEASE=NO`、`PROPOSED_RELEASE=UNSIGNED`、`D_TASK_CONTINUE=HOLD`。用户已验收原 1000 步并明确暂不继续；本交付不改变该决定。

## 已复核的 D checkpoint

从 D 的已归档分支提交 `145ee8050cf95a3a1988e868c6b935e10538f723` 只读取得归档。原始归档最初记录提交为 `11803909226538ab83847d1f06827a8292501139`。按实际文件重新计算：

- tar.gz：100,700 bytes，SHA-256 `b62a818e58a2cf423e4cb1b9ed0468ae73ad578a75c4b2cb9fb3f28baf75e7b1`；
- 外层 `archive-manifest.json`：SHA-256 `ab9f0dffa8ef145eac7792dfe1e14e07c1354df1c2152dafec97bad1dc3a57b1`；18 个 manifest 成员逐个大小与 SHA 相符；
- step-1000 `state.json`：SHA-256 `e08574678b3f8231c0a52ab25febb13c9d06e815b28149eb306d4b451c83ef6f`；数组 NPZ 和原子轨迹也与归档 manifest、checkpoint bundle SHA 一致；
- `source_manifest.json`：SHA-256 `9a694ea8f93655efa7a1c8695669ab2f17520c2d1b90f6df28115b2b1f9c2ad6`；
- NPZ 中为 64 个 Ag 原子，含 64 个 atom IDs、positions、momenta、velocities、masses、cell 与 PBC。ID 顺序与来源 manifest 完全一致；立方晶胞和三方向 PBC一致；全部状态数组为有限值；三条 packing/velocity/thermostat 随机流均为 PCG64。

完整机器结果见 `SOURCE_ARCHIVE_VERIFICATION.json`。此检查从 D 已提交归档字节重新验证身份，不会恢复或改写原 checkpoint。

## 执行器候选

`segment_guard.py` 是隔离代码，不修改原采样器。它把唯一允许边界写死为 global step 1000→2000 和追加 1000 步；每次积分前运行边界检查，step 2000 后不能取出 step 2001。执行前重新计算配置、代码、checkpoint tar、manifest 和本机势文件的 SHA，并核实任务、owner、checkpoint 和结构身份。旧 `--stage continue` 默认全局终点 120000 被视为禁止路径。

测试接口采用功能注入一步回调。每一步先在复制状态上推进并检查有限性；异常、超时、取消或归档失败即停止，不替换最后有效 checkpoint。完整 checkpoint 在新目录中 staging，写入内部 `MANIFEST.json` 与 `SHA256SUMS.txt`、逐文件 fsync 后原子发布；已有目录和文件不可覆盖。后续归档限制单文件 50 MiB，fresh-fetch/readback 必须逐层通过 SHA 才可报告持久化。

生产集成仍未完成：代码没有连接 ASE Langevin、没有真实的 ASE restart adapter，也没有在 D 的执行器中实现实际每步 stop hook。整套 gate 仅以脚本接口提案存在，不能在 D 上调用它来运行真实计算。将来集成需重新审查执行器，更新代码和配置 SHA，并单独签署。

## 合成测试

当前 19 项合成测试覆盖精确边界、step 2001 拒绝、PCG64 状态重现、身份/配置/执行代码/lock/势 SHA 的篡改拒绝、默认 120000 终点拒绝、超时/中断/NaN 停止、单步越过墙钟即丢弃该候选状态、atomic checkpoint 完整目录发布与禁止覆盖、50 MiB 门槛、内部归档 SHA、临时本地 bare Git push + fresh clone 读回，以及归档/远端失败时保留原 step-1000 checkpoint。测试仅用虚拟状态/虚拟一步函数；势文件哈希正向路径用合成 mock，不读取实际 EAM 数据。

临时本地 bare Git 通过不等于 GitHub 服务端 push 或实际任务分支恢复通过。后续执行前 D 必须分别在新 Codex 会话中对唯一归档做真实 GitHub 分支写入、fresh fetch 和内外 SHA 核验。

## 热化解释

D 原有 1000 步段仍是 0.5 ps 的 1600 K 初始热段，不是 1250 K 生产采样。原报告记录峰值 4266.30 K、最后 200 步平均 1838.35 K、负温度漂移、step 0→80 势能下降 31.67 eV、最近原子距离 2.20069 Å。这个状态支持存在明显初始瞬态且尚未显示数值不稳；它不支持“已热化”或“已成液”。64 原子有限尺寸下瞬时温度自然涨落应结合自相关修正窗口统计判断。

如果未来确需最小趋势验证，应另行批准从该确切 step-1000 checkpoint 运行到 global step 2000；仍保持原势、随机数状态、晶胞、温度、dt、摩擦和种子，最多追加 1000 步/0.5 ps；墙钟上限建议 600 s。step 2000 只交给 C 离线复审，不自动接着运行，不自动授予其他五条任务。

## 资源与停止条件

按上一次同任务实测，1000 步耗时 200.31 s、峰值 RSS 102,563,840 B、磁盘快照余量 28,792,524,800 B。原性能 JSON 的 `archive_bytes=172,708` 是本地输出文件树大小；外层压缩 tar.gz 实际 100,700 B。暂按最多 600 s 墙钟、1 GiB 可用内存余量、2,000,000 bytes 临时新增空间估算，但这些都不是 D 当前实时状态。获得下一阶段批准后须在 D 重做资源预检。

结构、势或原始归档 SHA 变化，配置/代码 hash 不一致，checkpoint不是 step1000，owner/task错位，stop 参数缺失，RAM/磁盘/inode不够，首个积分或写盘异常，出现 NaN/Inf、任意原子对距离低于旧执行门槛 1.8 Å、墙钟超限，或者任一远端上传/新克隆读回失败：立即停止。不得改势、改盒、重抽速度、重置 RNG、清理唯一有效副本、或用默认 120000 步继续。

## 交付状态

冻结候选 SHA-256：`segment_guard.py` `0da4ea9b02c948ed36e4e2f6a1a35edc8ba81dae537833d485a4cc033611084e`；`PROPOSED_RELEASE_1000_2000.json` `bc4c67a7227ce4006fb72ec63d57a08f3e4e7d907230b20f324deee64eb3a67e`；`RELEASE_LOCK_UNSIGNED.json` `bc3cae1e4f4f96f31aa07824c94441ad9b281fa2ff06ce0ff315c9875433212d`。所有交付文件 SHA 列在 `SHA256SUMS.txt`。

配置当前明确为 `PROPOSED_UNSIGNED_NO_MD`，release lock 为 `UNSIGNED_HOLD`；实际门禁使用本目录文件和真实归档执行后，在读取势文件路径前返回 `RELEASE_NOT_SIGNED_BY_C`。合成夹具中临时构造的 signed 状态仅用于测试 gate 的正向分支，绝非 C 签署。
