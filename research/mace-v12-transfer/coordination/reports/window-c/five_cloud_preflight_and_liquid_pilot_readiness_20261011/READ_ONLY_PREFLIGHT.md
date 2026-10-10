# owner补证清单

只在本机允许读取的环境内操作；不读A/B或其他主机的工作区。实时数据记录时间、repo HEAD、branch、命令及回执，脱敏后提交独立preflight报告并附SHA。报告无密钥、私人home路径或真实IP/machine-id。

| 项目 | 必需证据 |
| --- | --- |
| CPU | visible/nproc、affinity/cpuset、cpu.max和父cgroup限制、provider allocated CPU；最小quota而非5个可见vCPU |
| RAM | memory.max/current及父层、MemAvailable、effective headroom；E报available大于32GiB限额，不能当cgroup余量；H读数需说明来源，F/G补精确cgroup bytes |
| 磁盘 | 在最终work根挂载上df -B1/statvfs bytes、总量、inode、卷/用户配额、时间；五台现有值均低于30GB |
| 独立性 | 脱敏allocation/VM/node/volume标识及provider专用资源声明；不同namespace不足 |
| 持久性 | overlay backing volume和平台生命周期说明/备份恢复政策；不重启任何主机 |
| runtime | 文件路径、bash -n、内容审读、SHA、fresh shell的setup与MPI结果 |
| 软件 | 已有Python、GPAW26.7、ASE3.29、PAW1.2.1；gpaw info MPI/ScaLAPACK；OpenMPI5.0.7与GPAW链接ABI |
| setup | Ag/Ti/Si/C PBE.gz压缩文件SHA对照项目公开方法manifest，许可来源 |
| 任务 | 本机运行/排队ID与资源；不停止或重新分配现有任务 |

只读参考命令：os.cpu_count/sched_getaffinity、cat /proc/self/cgroup及可见cpu.max/memory.max/current、df -B1、mpirun --version、现有env的gpaw info、gpaw_data.datapath()。命令不足以证明父层/provider约束，owner还需补证。不要执行安装或H2 smoke。

四元素项目压缩setup SHA：
- Ag bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97
- Ti bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f
- Si 2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0
- C fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983

缺失项目为UNVERIFIED；不得安装/运行科学任务来填补空项。
