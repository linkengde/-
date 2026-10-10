# 独立持久目录和 runtime 方案

本文件为提案，未在D–H部署。WORK_ROOT_PLAN.csv逐台固定工作根、持久配置位置及branch。运行目录位于仓库外，不使用/tmp保存checkpoint；overlay不自动等于跨重启持久卷。

每台由owner核验 realpath、挂载来源/持久策略、volume配额/可用字节、写权限和脱敏provider allocation证据；不进入A/B目录验证。不能把不同路径或hostname当作物理资源独立证明。若/workspace是临时overlay，需在平台确认的持久volume上配置同等独立命名空间后再申请放行；本轮不挂载新卷或搬迁已有任务。

固定任务ID下按mesh_id分开inputs、SCF run、checkpoint和日志；同一host初批最多一个4rank SCF。checkpoint不得跨task/mesh复用。活动任务不搬迁、不停止、不重分配。sampler restart和DFT checkpoint分别归档。

持久runtime仅设置环境，不能自动安装、计算、后台启动或清理。F owner另提/workspace/.work-f/dft与/workspace/.work-f/sampling，尚未创建；也是可审候选。最终只登记一个DFT_WORK_ROOT，不并存活动根或搬迁已有任务；表中统一根是C提案，未远端应用。将以下设置写入本机已审查配置的提案需owner核验后应用：
- CLOUD_ID为d/e/f/g/h之一，DFT_WORK_ROOT为表中唯一值。
- GPAW_PY沿用现有安装的绝对interpreter，候选/workspace/.venvs/gpaw-mpi/bin/python，存在性先查。
- GPAW_SETUP_PATH由现有gpaw_data.datapath()唯一解析，四元素压缩setup SHA匹配；不追加未验证的setup路径覆盖优先级。
- 依据D实测build，GPAW_MPI_BACKEND=cgpaw；F/G/H须确认其现有build是否同样需要。不能假设不设置时裸Python成功。
- OMP_NUM_THREADS=1、OPENBLAS_NUM_THREADS=1、MKL_NUM_THREADS=1；MPI的PATH/LD_LIBRARY_PATH只沿用各自已核验ABI，不照搬其他主机library路径。
- runtime在/workspace/cloud-config/cloud-<id>/runtime.sh持久记录。现有D /workspace/onboarding/runtime.sh、E /workspace/.setup/e/runtime.sh、F /workspace/.setup/onboarding/runtime.sh、H /workspace/.setup/runtime.sh不直接覆盖；先审读/哈希，G可在确认activate.sh仅设置环境后引用它。

F/G/H检查重点：脚本存在性、bash -n、文件SHA、内容是否纯环境配置；fresh-shell setup解析；gpaw info的MPI/ScaLAPACK、实际扩展链接的libmpi及launcher ABI。先看文本再source未知脚本，避免隐含科学启动。已有版本完整不代表fresh shell可用，也不需要重装。现有安装脚本非check路径包含H2 DFT smoke，当前禁止运行。

如用户/owner允许非科学MPI sanity，可用已核验Python仅import gpaw.mpi并assert world.size=4，不构造calculator或调用能量/力；C本轮未远程执行。MPI失败报告，不自动oversubscribe、改backend/build、安装或绕过资源限制。

持久性由平台卷策略/备份与恢复约定证明，不通过停机/重启做测试。input/manifest/compact结果/校验可提交owner分支；full轨迹和checkpoint保留各自持久卷，禁止GitHub巨型checkpoint推送。任何清理或retention策略修改必须预先获准。
