# F 云端持久工作区与资源补证

采集时间：2026-10-11 01:41:05–01:44:03 UTC。OWNER=F，身份依据为当前会话的 F 指派、固定执行包内 `START_F.md`，以及用户补传文件与仓库文件 SHA-256 完全一致。该身份为 owner 自报，未由 C 远程核实。

执行包固定在 GitHub `linkengde/-` 提交 `2959fd3e2396a7319bd1f5c7f4fec27ce442abe0`。对该包运行 `sha256sum -c SHA256SUMS.txt`，27 项均为 `OK`；清单自身 SHA-256 为 `56522eaa220ac9a3c7ef87c9947f3a6a68ad4e35ffe699a7c4e399ada87866e1`。`START_F.md` 的 SHA-256 为 `415d0a59a8478bd26fe1919259fab25b55537cdd94a79b5c3f35cb0399ef7d14`，与补传附件一致。

| 项目 | 状态 |
| --- | --- |
| `RUNTIME_PASS` | YES（脚本语法、环境解析、版本、PAW 哈希及 MPI 动态链接静态核验通过；本轮未启动 MPI rank 测试） |
| `PERSISTENCE_VERIFIED` | NO |
| `SAMPLING_STORAGE_READY` | NO |
| `DFT_STORAGE_READY` | NO |
| `RESOURCE_READY` | NO |

当前没有可确认的持久卷。`/workspace` 是可写 overlay，未见平台配额或正常会话重启保留证明；`/tmp` 是 tmpfs。两个专属根目录都不存在，且没有被创建：`sampler_work_root` 与 `DFT_WORK_ROOT` 仍未配置。只有平台确认持久挂载后，才建议分别采用 `<verified-persistent-mount>/cloud-f/liquid-ag-six-pilot` 和 `<verified-persistent-mount>/cloud-f/ag-ti-si-c-b1`。C 配置中的 F 根目录、持久性和隔离放行字段仍为 `null/false`，本报告没有修改配置。

`preflight_owner.py` 对两个未创建的提案路径分别输出了原始 JSON；另对现存 `/workspace` 挂载点做了独立测量探针。提案路径报告保留 `work_root_exists=false`；`/workspace` 探针虽为 `work_root_exists=true`，仍明确记录 `persistence_verified=false`、`independent_resources_verified=false`，因此它只代表当前容器 overlay，不是已验证的 sampler/DFT 工作根。挂载、配额和三份原始 JSON 见相邻文件。

资源快照：CPU affinity/cpuset 可见 5 个逻辑 ID，但 cgroup `cpu.max=400000 100000`，有效配额为 4 CPU。cgroup 内存上限为 34,359,738,368 bytes，采样时已用 3,587,133,440 bytes，有效余量 30,772,604,928 bytes；系统级 `MemAvailable` 高于 cgroup 余量，故按较小的 cgroup 余量报告。F 当前可用工作盘为 28,873,228,288 bytes，可用 inode 为 2,036,751；provider/volume quota 未知。01:44:03 UTC 的后续原始挂载探针记录可用盘 28,873,207,808 bytes、可用 inode 2,036,747。容器进程清单未发现 GPAW、LAMMPS、MPI launcher 或本包采样/DFT runner 进程；没有可查询远端作业队列的接口。

环境复核：Python 3.12.14、GPAW 26.7.0、ASE 3.29.0、OpenMPI 5.0.7、gpaw-data 1.2.1、NumPy 2.5.3、SciPy 1.18.1、LibXC 5.2.3。`/workspace/.setup/onboarding/runtime.sh` 的 SHA-256 为 `299f595beb8233cc7b22601a01d1da7e63380ce4cc044d6cd67d8bf979a9ed00`，`bash -n` 通过；加载后 `GPAW_MPI_BACKEND=cgpaw`，`GPAW_SETUP_PATH=/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups`。`gpaw info` 显示 MPI 与 ScaLAPACK enabled；GPAW 扩展链接到 OpenMPI `libmpi.so.40` 和 `libscalapack-openmpi.so.2.2`，`ldd` 未报告缺失库。Ag/Ti/Si/C PBE setup 哈希见原始预检 JSON，四项均匹配项目基准。未进行计算或 MPI rank 测试。

**采样空间与 DFT checkpoint 分开估算。** F 自己负责的任务为 `C-F-AGL-1250K-rho100-seed02` 和 `C-F-AGL-1250K-rho103-seed01`。固定执行包中的 512 MiB/轨迹上限、32 MiB/checkpoint reserve、2 GiB minimum-free 条件目前仍为 proposal-only，且 `storage_release=false`、`warmup_release=false`、`continuation_release=false`。以 F 自己的两项任务计，保守相加预算为 1 GiB 轨迹归档 + 64 MiB checkpoint reserve + 2 GiB 空闲底线 = 3,288,334,336 bytes；当前 overlay 数值上足够。六项分布在 D/E/F 三台时，每台仍需单独至少该容量；三台容量不可池化，D/E 工作区证据缺失，所以不能宣告六项整体就绪。F 的 `SAMPLING_STORAGE_READY` 仍是 NO，因为目录、持久性、配额和独立性均未验证。

GPAW 写盘另按 C 的工程估算参考：Gamma 与 2×2×2 阶段累计需求示例为 16.773–23.402 GB；若扩到 3×3×3，示例为 31.032–50.422 GB。该估算假设未验证的密度、band 数和 checkpoint overhead，不是实测或容量保证。当前 28.873 GB overlay 对较密网格可能不足，且没有 checkpoint retention ledger；因此 `DFT_STORAGE_READY=NO`。不删除或搬动任何现有 checkpoint。

当前执行包 `C_RELEASE_CONFIG.json` 的 SHA-256 为 `57dc948d45daa3c51f7b3e5491a19e4f407546806531c3afae14a221b054ffc0`；F 的 `sampler_work_root`、`dft_work_root` 与 owner evidence SHA 仍为空，host persistence/isolation 仍为 false。没有修改任何放行字段。F 两个 task ID 保持原分配，仍等待 C 收齐证据、签署 host storage release 并发布暖机配置 SHA；在此之前保持 HOLD。没有安装软件、生成结构/轨迹、启动采样/MD/DFT/训练、访问 A/B 文件，或修改 C 主配置。

报告文件按本目录 `SHA256SUMS.txt` 校验。代码分支为 `cloud-f-b1-dft`，以指定提交为基点。
