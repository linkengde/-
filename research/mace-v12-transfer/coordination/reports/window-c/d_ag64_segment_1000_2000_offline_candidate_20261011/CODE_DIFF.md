# 与上一版 `segment_gate_proposal.py` 的差异

此前原型只检查若干由调用方传入的哈希字符串和身份字段，没有重算本地文件字节，没有读取归档成员，没有实际 stop-loop，也未接入生产执行器。

本目录新增的 `segment_guard.py`（SHA-256 `0da4ea9b02c948ed36e4e2f6a1a35edc8ba81dae537833d485a4cc033611084e`）：

- 对真实 step-1000 tar、外部 archive manifest、18 个归档成员、state JSON、source manifest、NPZ 的 64 Ag 原子数组、原子 ID 顺序、晶胞/PBC、PCG64 三条 RNG 流逐字节重算 SHA/身份；读取 D 归档为只读。
- 要求实际配置和执行代码路径，现场重算其 SHA，再与独立 release lock 中的 pin 比对；实际本机势文件亦要求重算指定 SHA。
- 将起点、终点和新增步数钉死为 1000、2000、1000；循环在每一步前拒绝越过终点，因此不能走到 2001，也不会触发旧 `--stage continue` 默认终点 120000。
- 只通过功能注入接口接受一步函数；每一步复制状态，有限性检查通过后才替换已验证状态。超时、取消、异常或 checkpoint 保存失败均终止。
- 通过临时目录、manifest、内部 SHA 清单与单次目录 rename 原子发布完整 checkpoint；既有路径拒绝覆盖。输出 tar 的每个成员限制为 50 MiB，并排除势文件；本地 bare Git 的 fresh clone/readback 只在测试中模拟。

仍未实现、也不能据此运行科学计算：没有将该协议接入 `run_sampling_codex.py`/ASE Langevin；没有真实 checkpoint 到 ASE 原子状态的生产恢复适配器；没有执行真实积分、D 端实时资源预检、GitHub 推送或远端科学结果归档。将来集成这些功能必须单独改代码、补测试、重新计算完整 SHA、C 审查并签署新配置。

新配置保持 `PROPOSED_UNSIGNED_NO_MD`，锁文件保持 `UNSIGNED_HOLD`。测试里的 `SIGNED_BY_C` 只是临时合成夹具，用来检验正向控制流，不构成签署或授权。
