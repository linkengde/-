# C：Ti₃SiC₂ (0001) 两终止面 DFT 输入审核

日期：2026-10-17（Asia/Shanghai）。E 输入包提交：`cbc14c97fe5d962dc808a4cf48bca3e6b5588ea0`。对应 C 任务交接：`1b641b31e750731509a34f033eb5eeee47c029d7`。

## 裁定

INPUT_STRUCTURE_REVIEW = APPROVED
SCIENTIFIC_METHOD_PROPOSAL = APPROVED_FOR_FIXED_GEOMETRY_SINGLE_POINTS
ROLE = FOUR_TRAIN_CANDIDATES_WITHIN_EXISTING_BULK_PARENT_LINEAGE
DEV = HOLD_NO_VERIFIED_INDEPENDENT_SOURCE
DFT_EXECUTION_READY = NO
execution_authorized = false

本裁定批准四个结构作为后续标签的 TRAIN 候选，不等于 DFT 执行授权或标签已产生。E 所述“用户已授权补充计算”在此具体落到四个已冻结表面输入；但资源/执行门尚未通过，故当前不能运行。

## 结构与来源

从 E 固定提交读取包内文件，14 项 SHA 清单均通过；其离线几何验证器报告 `PASS_GEOMETRY_AND_HASH_CHECKS`。C 用 ASE 3.29.0 独立读取四个文件，确认各100原子、PBC=TTT、晶胞向量长度6.115/6.115/53.247 Å、面内夹角120°、表面跨度35.247 Å、真空18 Å；Ti终止组成Ti52Si16C32，Si终止组成Ti48Si20C32。各参考-响应对只移动清单所列两个同种表面位点，法向位移为−0.02/+0.02 Å，原子顺序和晶胞不变。最近邻2.0882 Å，无短接触。输入 SHA 与 manifest 一致，坐标没有标注DFT结果。

来源为 COD 9009647 / Kisi (1998) 实验晶体，源 CIF SHA `6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3`。这是已有48原子体相来源的同一晶体学根源；四个候选只能留在TRAIN侧，不能被称为独立 DEV/TEST。两终止面不是来源独立的泛化家族。表面由实验坐标对称切割并重复，未弛豫；固定几何力可以为所列欠配位构型提供训练信号，但不据此声明稳定表面、表面能或平衡重构。

对于Si slab报告的镜像误差约3.06×10⁻⁵ Å，属于来源 CIF 坐标精度内的微小残差；响应两端的位移在当前原子坐标中仍严格为±0.02 Å。该残差不影响几何身份或输入资格，不把它包装成更高对称精度。

## 参数审核

批准固定几何 single-point 提案：PBE/PAW，PW500 eV，Fermi–Dirac0.10 eV，SCF energy/density/eigenstates阈值各1e-5，maxiter160，成功P0的MSR1参数，GPAW默认Davidson niter2，4 MPI / OMP1。结构不松弛、不改胞。surface应力仅保存诊断，真空胞的3D体积归一化应力不作为体相弹性标签。

MP `(6,6,1)` 标准整数网格的坐标定义清楚：同一面内晶格下，其6×6面内点与体相 `(6,6,2)` 网格投影一致，z仅为Γ。这个比较证明网格对应关系，不证明新 slab 的网格收敛。作为这四项力训练候选的初始网格可接受；结果只能按其实际有限网格标记，不能声称应力或表面能已收敛。

实际method/role执行锁必须以新结构SHA重写并逐项绑定，不能直接复用P0锁SHA。E包引用的旧锁和runner仅证明实现来源，不授权或冻结新任务。

## 角色与数据限制

下列四个未来标签的拟议用途获准为TRAIN，且只在其固定结构产生通过验收的DFT标签后生效：

- `TSC100_TI0001_COD9009647_2x2x2_REF` — `b1948a44b89efe17ee7c2891b984167982067feeffa774b7f6fa77d9e9158d86`
- `TSC100_TI0001_COD9009647_2x2x2_RESP_Z020` — `afa99ddae62e861bf641d3e7f2d908c78dfb0d466e762a25ebaeec64b0a91078`
- `TSC100_SI0001_COD9009647_2x2x2_REF` — `ac8b804516424215694c284d3bcc36d23da590f80390f47a0805c87cedb764c4`
- `TSC100_SI0001_COD9009647_2x2x2_RESP_Z020` — `f78560dfe891ea86d8a7da5843216050430f9aef8303c463231f31c19ca10da0`

每个参考/响应对须作为同一子族处理；与bulk_parent_48共享COD 9009647上游来源，不能拆作DEV/TEST。能量与应力保留诊断，训练只允许逐原子力目标。当前六个已批准体相标签仍是独立的有限探索训练映射；本文件不改写任何旧角色档案。

E提出的两个DEV候选没有结构/来源，维持HOLD。不得拿这两个表面家族冒充独立DEV。封存TEST仍由独立evaluator负责；C没有查看或枚举它。

## 为什么 DFT 仍 HOLD

E资源记录为4核、32 GiB内存、约16.66 GB可用工作盘。内存估算约6.18 GB波函数+1.32 GB投影子数组，不含FFT、对角化、密度混合及Python峰值，未对100原子slab实测。单点时间粗估约3.7小时；外部walltime尚未查明。10 GiB临时scratch预算未经实测，若占满只剩约6.66 GB余量。无经过验证的slab中断恢复/checkpoint流程。新任务输入、执行锁、runner和结果归档方案尚未冻结。

运行前的最小放行条件：E在其窗口确认单任务5小时作业窗口可用；实测空闲内存/磁盘并保留至少10 GiB计算scratch和6 GiB余量；一次非科学性runner/原子归档验证或明确“不checkpoint、失败留日志、不得自动重试”的停止路径；新锁绑定四个输入SHA、PAW SHA、MP坐标约定、runner SHA、资源与归档策略。四点严格串行，一次仅一个GPAW进程，RSS达到28 GiB或可用盘低于6 GiB即停止；墙钟达到5小时立即安全中止并归档失败日志。成功结果每项原子化归档、fresh-fetch并核SHA后再启动下一项。不得清理历史科研文件腾空间。

通过以上资源与执行门后，可按这四项既定输入执行，不需要另行申请扩大任务范围。当前没有启动DFT、拟合、MD或TTM。
