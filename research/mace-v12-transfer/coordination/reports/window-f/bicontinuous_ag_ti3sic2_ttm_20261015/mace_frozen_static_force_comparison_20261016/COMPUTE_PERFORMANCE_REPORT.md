# 全模型 MACE 静态推理：性能与完整性

| 指标 | MACE-MP | MACE V17 |
|---|---:|---:|
| 完整 tile 数 | 3997 | 3997 |
| 中心原子数 | 39486 | 39486 |
| 累计 tile 力推理时间 | 27.98 h | 19.58 h |
| 平均 tile 推理时间 | 25.20 s | 17.63 s |
| tile 推理中位数 | 25.83 s | 18.00 s |
| 单 tile 最长推理 | 49.61 s | 31.35 s |
| 峰值 RSS | 29.33 GiB | 29.04 GiB |
| 最大裁剪原子数 | 3262 | 3262 |
| 最大裁剪有向边数 | 190278 | 190278 |

两模型使用 float64、CPU、每模型串行推理；PyTorch 设 5 个线程。每个核心块按 Y/Z 周期最小镜像展开、X 非周期，并带 18.1 Å halo；MACE cutoff 为 6.0 Å。环境上限 32 GiB，实测峰值 RSS MP/V17 为 29.33/29.04 GiB。

累计 tile 推理时间为 3,997 个 tile 的 energy+force 计时相加，不含模型载入、邻居图构建、磁盘写入和进程中断等待；MP/V17 分别 27.98/19.58 小时，总计 47.56 小时。环境中断后均从已保存 tile 断点续跑；最终覆盖按每块预期 atom ID 核验，未见 OOM 记录。

源结构 SHA-256 e2b62a68c46e23037370a6ca20951aa0774edea2d728d3c3aadbcba4b3755674；MP/V17 权重 SHA-256 2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012 / 559b12c024fe8e2155e26ffc7ca48cf9fa1c2682ea37ff9c5367af14a5cc43d3。输出 MP/V17 全原子力 CSV SHA-256 d709ef07987a8706b733aa34a71934151fb3048a6168d6029b7041a96b4791e6 / 27f3bd8bf8a973d6aa631b73454b89bb93210f5b021f9dc831f202ce9aeecb29.

全量逐原子力及力差见 MACE_ALL_ATOM_FORCE_COMPARISON.csv；高力原子见 MACE_HIGH_FORCE_ATOMS.csv；区域和阈值敏感性见 MACE_FORCE_REGION_METRICS.csv。没有启动 DFT、训练、松弛、300 K 测试或 TTM。
