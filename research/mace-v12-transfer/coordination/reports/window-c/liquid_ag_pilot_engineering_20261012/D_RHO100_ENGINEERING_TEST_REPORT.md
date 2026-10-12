# D：rho100 Ag64 首段工程预检

状态：`ENGINEERING_PREFLIGHT=PASS`；C 的科学锁保持 `execution_authorized=false`，正式任务仍为 HOLD。预检报告、程序和测试均绑定在本目录；C 固定输入为提交 `287b60af424f7bcc3a7b3598e524f5b543884ec9` 的 `LIQUID_AG_PILOT_LOCK.json`，SHA-256 `e1adaf12d542184c709d2eb2a5363cbfe9578c168ebc5e7ad260ee39a4faf09b`。C 提交内的 `SHA256SUMS.txt` 对该锁文件的条目也匹配。

## 预检结果

- 本地 `Ag.eam.alloy` 为 780452 字节，SHA-256 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`。使用 C 固定版本的 `inspect_eam.py`（SHA-256 `d3a437582dcd5bad7e361c6c81f64ac1191161170105b43133ffdd4065526a77`）按 ASE 3.29.0 解析为单元素 Ag 的 setfl/alloy；ASE EAM 计算给出有限的 eV 能量、eV/Å 力和 eV/Å³ 应力。两原子力的中心差分误差为 `3.65e-11 eV/Å`。
- 锁定 cutoff 为 5.995011 Å，盒长一半为 5.358728 Å。半盒长边界测试保留了不同整数周期镜像偏移 `(-1,0,0)`、`(0,0,0)` 和 `(1,0,0)`；没有把镜像邻居折叠成单一原子索引。ASE EAM 镜像构型的能量、力均有限。
- D 的顺序拒绝初始化与 C 固定提交里的 `pack_random` 逐位比对通过：位置数组及三个 PCG64 子流状态完全相同，均在 2598 次候选点后放入 64 个原子。临时测试构型的最小周期距离为 2.204983 Å，原子 ID 为唯一的 1–64。该临时构型仅在 `/tmp` 使用，测试结束后已清理，没有加入本报告包。
- 温度初始化使用锁定的 ASE MaxwellBoltzmannDistribution、`force_temp=True`、1600 K、velocity 子流和 FixCom。ASE 计算自由度为 189，约束感知温度为 1600 K；Langevin 使用锁定的 0.5 fs、0.01 fs⁻¹、`fixcm=False`。温度和能量日志字段均为有限数值。
- 在两个相同的临时构型上比较不间断 6 步与 step 3 接收 SIGTERM、保存后恢复到 step 6：最终坐标、动量、三个 RNG 状态和 0–6 的全部轨迹帧逐项完全相同。恢复前核验了 checkpoint、trajectory、thermal log、energy log 的 SHA 和文件长度。SIGTERM 在已完成步边界停止；没有生成正式轨迹或正式 checkpoint。
- step 1001 和墙钟 3601 秒的守卫测试均按锁定上限拒绝；边界 1000 步和 3600 秒通过。runner 的正式路径只有在 C 单独签发并提供匹配 SHA 的授权 JSON 后才会继续，当前 C 锁不会触发正式结构或 MD。
- 只用重复静态测试帧（0 个 MD 步）测量 1000 步锁定输出频率对应的归档体积估算为 59,608 字节，最大单文件远低于 50 MiB，低于 512 MiB bundle 上限。当前 `/workspace` 可用 28,785,012,736 字节、inode 余量 2,037,659；cgroup 配额为 4 CPU、内存上限 34,359,738,368 字节，预检进程峰值 RSS 约 107 MB。2 GiB 空间门槛和 32 MiB checkpoint 预留均满足。

ASE 的 C 参考解析器将径向尾部状态标为 `C_TAIL_REVIEW_REQUIRED`，这是解析器不签发策略许可的状态；同一份已验 SHA 的 C 锁明确采用仅适用于该原始文件哈希的“一网格 bin 近零尾部”策略。解析输出确认尾部范围在一个网格间隔内、数值约为 (10^{-13})，D 没有更改势函数或自行替代 C 的策略判断。势函数文件本身及论文 PDF 不会上传到公开分支。

## 绑定摘要

- runner（`pilot_engine.py` 与 `run_pilot.py` 组合 SHA-256）：`af9c3a0b223b5506e6127af73d02387eec4000580089ab2f7f584fbca793a9f8`
- generator（`pilot_engine.py` SHA-256）：`ce3b4329b6c650398a422fbe0ef94b949ed83c073ff3c59d57e6bdcfc40e47a8`
- environment（Python/GPAW/ASE/NumPy/SciPy/OpenMPI、runtime.sh 与平台信息的规范化 SHA-256）：`f3f6d27f7e8b387bb5cc7b91048dba7fa26280f9ed54aa1572c73879946e867f`
- runtime：Python 3.12.14、GPAW 26.7.0、ASE 3.29.0、NumPy 2.5.3、SciPy 1.18.1、OpenMPI 5.0.7；`GPAW_SETUP_PATH` 与 `gpaw_data.datapath()` 一致。
- 六步临时工程段的测得时间约 1.57 秒；按该测试的逐步 checkpoint 和文件 I/O 线性外推，1000 步约 262 秒。这只是资源估算，不是正式运行时间承诺。

## HOLD 与提交范围

正式初始结构未生成；正式 MD、DFT、MACE 训练均未启动，正式步数为 0。step8000 的 `ARCHIVED_HOLD`、其他五条任务、科学参数和数据角色均未修改。

上传内容只包括本工程预检代码、C 的元数据检查脚本副本、测试和报告，不包含 EAM 势函数、PDF、正式结构、轨迹或 checkpoint。D 分支推送后会在全新临时 clone 中逐文件校验 SHA；校验结果及提交 SHA 由本次交接信息提供。

下一步由 C 审核这些固定 SHA，并单独签发匹配本任务、锁文件、runner、environment、1000 步和 3600 秒上限的 `execution_authorized=true` 授权。未经该授权，本 runner 保持 HOLD。
