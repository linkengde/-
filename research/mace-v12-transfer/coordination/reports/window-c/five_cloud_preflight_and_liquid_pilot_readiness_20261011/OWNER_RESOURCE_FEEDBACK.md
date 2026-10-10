# Owner 实时反馈的来源记录

本文件将本轮用户转交的五台实时反馈规范化保存。原始命令回执/文件SHA/测量时间未提供，未形成五个host自身Git报告；不等同C远程检查。每位owner对“其他主机UNVERIFIED”的声明只限其可见范围，C汇总时采用各主机自己的反馈，不用一台代填其他台。

## D

branch cloud-d-b1-dft，repo snapshot f9544cd00046e7fd4be75fec4d6fe13e50e8cdab，owner报clean。cpu.max400000/100000=4，nproc5。memory.limit34359738368、cgroup余量30804070400 bytes；workspace total33770192896、available28797075456 bytes、overlay。
/workspace/onboarding/runtime.sh语法PASS、GPAW Python存在，setup解析到/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups。不设MPI backend的裸Python探针失败，GPAW_MPI_BACKEND=cgpaw后4rank成功；脚本未source/修改。owner RESOURCE_READY=NO。

## E

cpu.max400000/100000=4。报available34777948160、limit34359738368 bytes；由于available超过limit，不能当成cgroup剩余额度，须核对memory.current/父层与系统MemAvailable。workspace available28800167936 bytes。
/workspace/.setup/e/runtime.sh语法和setup检查PASS，独立持久work root BLOCKED。未附repo报告SHA。

## F

owner当前work/bcb51e5，clean；交接固定bb8不属于其所读origin/main684e历史。4vCPU quota、可见5logical IDs、physical allocation UNVERIFIED。cgroup limit32GiB、current约3.36GiB、headroom约28.64GiB；系统MemAvailable33931336kB为另一个口径。workspace available28874539008 bytes。
Debian13.6、Python3.12.14、GPAW26.7、ASE3.29、OpenMPI5.0.7、gpaw-data1.2.1、numpy2.5.3、scipy1.18.1、LibXC5.2.3。MPI/ScaLAPACK enabled，动态库解析无missing，历史4rank验证PASS；本轮未运行计算。
/workspace/.setup/onboarding/runtime.sh可加载，setup路径和四元素hash与项目基准匹配；package license GPL-3.0-or-later为owner报告。DFT_WORK_ROOT unset，owner提/workspace/.work-f/dft、/workspace/.work-f/sampling，尚未创建。
LAMMPS未安装。当前准备路线优先ASE/MACE结构采样候选与GPAW，LAMMPS不是本阶段依赖；不因其缺失否定DFT核心栈或发起安装。owner总体ENV_READY=NO包含更广范围，保留其原述但不作为核心栈重装依据。

## G

Debian13，4CPU quota400000/100000，allowed CPU集合5，纠正此前5核称谓。limit32GiB，headroom约28.9GiB、exact bytes未给。workspace available29018927104 bytes。
runtime.sh未找到；/workspace/.setup/activate.sh此前setup与gpaw_data.datapath一致，fresh shell回执未给。work/bcb51e5，clean。不把交接bb8当成本机实测报告。

## H

owner当前work与本地origin/main均bcb51e5。4.00CPU等效配额、guest可见5vCPU；原cpu.max/affinity回执未给。reported available34204598272 bytes，需系统/cgroup口径与limit/current补证。workspace available28896153600 bytes，overlay跨重启持久性UNVERIFIED。
/workspace/.setup/runtime.sh，setup指向gpaw_data.datapath且目录存在，mpiexec/mpicc OpenMPI5.0.7，GPAW links libmpi.so.40/ScaLAPACK、无缺失库。file SHA和完整文本未给。

五台软件/四setup总体匹配声明与上述逐机回执共同记录为OWNER/USER_REPORTED；provider独立性、持久性、quota与实测checkpoint大小仍需补证。没有任何结构或DFT执行授权。
