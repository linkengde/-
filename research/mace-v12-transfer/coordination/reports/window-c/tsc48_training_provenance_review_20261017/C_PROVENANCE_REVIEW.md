# C：六点 P0 来源补证复核

日期：2026-10-17（Asia/Shanghai）。E 固定提交：fa699e1340e8777a7095bb8c639a7247c28d1956；前序标签审核：b3eef435e372b833053b369b4d401165be28f745。

## 最终裁定

PROVENANCE_REVIEW = ACCEPTED_WITH_LIMITATIONS

TRAIN_USE_APPLICATION_READY = YES

DATA_ROLE = DIAGNOSTIC_ONLY

TRAIN_USE = HOLD；TRAIN=[]；DEV=[]；TEST=[]。

本次来源补证已足以进入明确范围、声明数值局限的父族训练用途申请，不要求继续追索历史上未采集的逐作业 PAW 哈希，也不要求仅为补 launcher 退出码重算。提交申请并不等于用途获批；本轮不分配训练角色、不拟合、不运行计算。

## 核验结果

真实 GitHub fetch 固定提交；按仓库根相对路径读取补证清单，16 项 SHA 全部通过。六份 COMMAND 的输入 SHA 与现有输入逐项相符，锁及 runner 原始文件字节哈希与 COMMAND、批次预检、补证 manifest 一致。执行锁 SHA=81d9ce869d153b3bde909f3d8259a3b89a437a8150158dcd32211e9a911afd0f；runner SHA=812110a1d63c8b374f5ace669fec00aac26b04ae125dccbe102ad197382cc9e4。

代码静态审查确认实际 GPAW 参数来自 scientific_method 中的 kpoints.size=[6,6,2]，没有执行任何6×6×4说明字符串；固定几何，500 eV、PBE、FD0.10 eV、MSR1、SCF三阈值1e-5，输出 native energy、force-consistent free_energy、完整 forces/stress。未导入或执行该科学 runner；复核脚本只读取 Git blob 和已有数据。

旧 MESH_SELECTION_DECISION.md 及锁内 mesh_selection_evidence.limitation 的6×6×4文字属于过期非执行说明。依据有效锁字段、runner参数路径和前次已审核的六份原始日志，正式裁定实际六点均为标准 MP6×6×2。保留历史文件与原SHA，在新增记录中解释，不伪装为6×6×4数据。命令本身只绑定锁，网格由锁和实际日志确定。

## PAW证据强度

Ti/Si/C 的批次预检、E_RUNTIME_PREFLIGHT、候选来源 manifest 及本次PAW哈希声明一致，与项目基准一致。六项日志内部PAW ID的一致性已在前次C审核确认。本次C验证的是这些归档记录之间的一致性；当前E本机文件哈希由E声明，C未远程读取其会话文件。

runner只检查GPAW_SETUP_PATH存在，不逐文件哈希；锁没有PAW SHA字段；批次预检的execution_scope只列VOL_m02。故这不构成六作业或每rank的历史现场哈希证明。作为现有诊断及未来有限探索用途申请的来源证据，可接受这项已知限制；不得写成逐作业已验证。后续新计算应在启动前绑定PAW SHA，本轮不追溯改写旧锁。

## 网格证据及适用范围

独立读取母体6×6×2和6×6×4原始extxyz，校验其SHA并比较元素顺序、晶胞及坐标；使用已有force数组与energy及summary中的free_energy复算：

- native差：+0.3183363911 meV/atom；
- free-energy差：+0.1578549395 meV/atom；
- 逐原子力向量RMS差：0.0028842303 eV/Å；
- 最大逐原子力差：0.0048810019 eV/Å。

符号均为6×6×4减6×6×2，与补证相符。该对照只验证母体上增加kz的效应，不是相对于无限密网格的总误差上界，也不证明六个扰动几何或应力已收敛。沿面内方向网格加密的证据需另行按其实际比较范围解释，不能用这个kz对照单独认证。

现有六点足以开展已获准范围内的冻结势诊断；如申请有限探索拟合，需明确接受母体网格精度可近似转移到小扰动结构的未验证假设。不得据此宣布弹性常数、应力约束或生产势精度通过。若拟合目标要求严格收敛的应力/弹性或更高精度标签，应另立最小收敛核验提案，不能自动增加DFT任务。

## 训练用途申请建议

建议先提交六点的有限探索性、力重点训练用途申请，暂不启用应力loss，也不据此标定弹性常数；保存应力供诊断。能量保留双口径；若未来加入能量loss，建议与本批force/stress匹配的有限展宽free_energy，明确是0.10 eV电子展宽自由能，不能与其他数据的extrapolated energy混用，也不能把电子展宽直接称为物理离子温度。此为用途申请建议，未授权实现。

bulk_parent_48整个父族必须作为同一来源单元管理；批准六点不自动批准旧五点或H候选。现有旧数值控制档案保留原角色，任何额外训练复用须在用途映射中明确逐项批准。同族所有后代不得用于独立DEV/TEST。独立DEV仍需另选来源；TEST由独立evaluator管理，本轮未读取或枚举封存数据。

附 TRAIN_USE_APPLICATION_DRAFT.md 可直接用于下一次用途审批。没有创建训练数据集、复制标签或启动任何科学计算；拟合集仍为空。来源审核闭合后不再要求E反复发布同类补证。
