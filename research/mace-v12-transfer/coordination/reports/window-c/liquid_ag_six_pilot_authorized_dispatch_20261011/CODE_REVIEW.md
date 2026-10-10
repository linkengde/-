# C执行包离线验证与独立复核

PACKAGE_VERIFICATION.json记录30项根C实际检查：七脚本AST／标准库help、JSON/schema、六任务身份／seed／owner、density-mass-cell、source/parameterSHA、config绑定、主机fingerprint算法一致、默认关闭门不启动科学作业／不产候选、错误配置SHA拒绝。

C另行独立评审发现并修复：MSD高R²仍可能是弹道（加入晚时指数和三窗口斜率一致性）；density gate遗漏；symlink路径隔离；checkpoint文件与目录fsync顺序；累计wall／暂停续算许可复用；官方setfl末端单dr尾段显式policy；资源budget重签与科学身份过滤接口；MPI初始化／非零退出暂停后续。

复核范围与证据：

- 采样／QA：density=false拒绝、owner路径逃逸／symlink拒绝、原子JSON持久化纯文件测试、resource-only预算重签与dt科学变化区分。线性MSD通过，弹道alpha=2／ratio≈2.98及平台MSD拒绝，未冻结阈值HOLD。仅合成数值数组，未生成原子构型。
- 重启：静态检查精确arrays／RNG、完整checkpoint文件和目录fsync后发布最新state；每stage累计wall、unclean停机保守收费、resume绑定当前checkpointSHA。无真实MD重启测试，首次授权暖机后仍须验证真实restart身份／性能。
- DFT：独立核对GPAW26.7源码hook/collective checkpoint/API；纯mock证实资源失败／子作业非零退出暂停后续，False模板不写ownerroot。无MPI或SCF实跑，首次数值control仍需真实验证。
- 实际源：独立逐字节核对4参数SHA／大小与setfl完整有限表；纯SciPy尾段11点重算吻合，未计算原子能量／力。近零不称为数学严格零，source允许不自动放行存储或科学任务。

当前未发现阻止此离线执行包冻结的问题；这不是DFT／MD验收。SOURCE_PASS=YES，STORAGE／WARMUP／CONTINUATION／STRUCTURE／DFT均尚未放行。代码后续变更须更新SHA并重新复核对应部分。完整轨迹／GPAWcheckpoint没有在C运行。
