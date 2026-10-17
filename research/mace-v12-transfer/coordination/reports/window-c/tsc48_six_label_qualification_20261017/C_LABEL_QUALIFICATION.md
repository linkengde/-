# C：六点 P0 DFT 标签资格审核

日期：2026-10-17（Asia/Shanghai）。审核 E 固定提交 `f9f2e88e34b4c1e18e6a42ed16fd7c949f7c22ef`；前序父族用途审查为 C commit `de24ee62e38b33bae87db0cfb4c6198e84cd4e95`。

## 结论

LABEL_QUALIFICATION = PASS_WITH_LIMITATIONS_FOR_DIAGNOSTICS

DATA_ROLE = DIAGNOSTIC_ONLY

TRAIN = []；DEV = []；TEST = []。

TRAIN_USE = HOLD_PENDING_NUMERICAL_AND_PROVENANCE_REVIEW_AND_EXPLICIT_USE_APPROVAL

START_FIT = NO；NEW_CALCULATION_AUTHORIZED = NO。

六个已有结果通过基础标签完整性与收敛资格审核，可作为已收敛的诊断参考归档。不需要仅为补齐 launcher 退出码而重算。这里的 PASS 不代表新增 k 点精度或应力精度已经得到验证，不授予训练用途或任何新计算权限。

## 独立读取及核验

从真实 GitHub fetch 固定 E 提交，将该结果目录导出到 C 临时区。审核程序仅用 ASE 读取现有文件及 numpy 比较数组，没有构建计算器或运行 SCF、run 0、积分、拟合。

六个结果目录及审计目录内 51 条 SHA256SUMS 记录均通过；发布 manifest 中的输入、labelled.extxyz、RESULT.json、RUN_STATUS.json 和每项清单 SHA 均逐项核对。六个输入 SHA 与先前 C v1 manifest 相符。

每项均为 Ti24Si8C16、48 原子、三维周期；输入与输出的元素顺序、坐标、晶胞、PBC、local_new_id 和 spacegroup_kinds 完全一致。所有结果包含有限的48×3力和3×3对称应力。

原始 gpaw.log 各自确认 Converged in 16/17/16/17/24/24 steps。RESULT、发布 manifest、labelled.extxyz 的能量和力一致；原始日志的能量、力及应力在打印精度内一致。应力 eV/Å³、GPa 换算与 ASE Voigt 次序 [xx,yy,zz,yz,xz,xy] 已交叉核验。全精度 force 数据以 RESULT.json 为准。

标签中的 energy 对应 GPAW 日志的 Extrapolated 能量；free_energy 对应有限展宽的 Free energy，两种口径分别保存。用于能量—力联合拟合时必须先冻结采用的口径；不能静默把两者混用。

## 缺失 launcher 退出码

一项 exit_code=0，另外五项 supervisor 退出码缺失，仍记录为未知；不补写为0，也不宣称全部进程正常退出。原始 SCF 收敛行、最终能量、完整 forces/stress 和相互一致的静态结果文件足以支持诊断标签的科学完成状态。该限制可以接受用于现阶段诊断，不能代替程序退出/运行生命周期的完整溯源。

## 实际科学设置及训练前限制

原始日志确认 PBE、PW500 eV、Fermi–Dirac 0.10 eV、SCF energy/density/eigenstates 各1e-5、maxiter160、标准 Monkhorst–Pack 6×6×2、MSR1。RESULT 记录 GPAW26.7.0、ASE3.29.0、默认 Davidson niter2；COMMAND 记录4 MPI ranks、OMP1。六项 PAW 日志内部 ID 一致，见 v2 manifest。

这与最早的10×10×4、Pulay/PPCG提案存在差异。混合器/求解路径不同本身不使已收敛标签无效，但不能声称沿用原提案执行；特别是6×6×2对这六种变形的能量、力、应力收敛尚未在本次发布证据中证明。旧母体8×8×4与10×10×4的通过结果不能替代6×6×2证据。未来如需统一精度，先审核现有网格对照并量化差异；本轮不要求或授权重算。

COMMAND 绑定执行锁 SHA `81d9ce869d153b3bde909f3d8259a3b89a437a8150158dcd32211e9a911afd0f` 及 runner SHA，但相应 EXECUTION_LOCK_K6X6X2_MSR1.json、run_single_point_msr1.py 未出现在本次固定提交所声明的路径。因此已核对的是这些字段的一致性，不是锁/runner原始字节的独立复算。现有材料还需补齐实际运行 PAW 文件 SHA 与锁的关联，不能仅用日志路径代替文件身份。最小补证是发布现有锁、runner及PAW哈希元数据和已有收敛对照，不需要重跑结果。

labelled.extxyz 仍带旧 dft_label_status=PENDING_NO_DFT_RUN 字段，与已有结果不符。这是输入元数据继承问题；本轮保留文件字节，由新版用途 manifest 明确解释实际标签状态，不修改原件、不以该旧字段否认日志和标签。

## 父族与后续动作

所有六点仍属于 bulk_parent_48，同族 baseline、Ti/C位移及其他后代不能作为独立 DEV/TEST。新版 FAMILY_USAGE_MANIFEST_v2.json 仅更新标签证据状态，不变更任何历史角色。

下一步：E 补充现有执行锁、runner和PAW绑定证据；C 审核6×6×2现有数值精度证据及拟合能量/力口径；然后提出明确的父族训练用途申请。可在另行获准的诊断任务中以这些标签比较冻结 TSCBOP，结果仍为诊断，不作为独立验证。训练资格与诊断计算权限分别处理。

E6剩余8点仍是特定覆盖方案下的条件预算，不因本次诊断标签通过自动获批。独立 DEV/TEST 与clean surface的来源资格仍须分别审核，TEST由独立evaluator管理。当前没有新增DFT任务，没有读取封存结构或标签。

## 交付及复查

SIX_POINT_LABEL_REVIEW.csv 给出逐点数值和结论；FAMILY_USAGE_MANIFEST_v2.json 保存结构/结果 SHA、缺失退出码、PAW日志ID及用途状态；SOURCE_SHA256.json 保存读取的结果文件哈希。audit_published_labels.py 可在导出的 E 结果目录上复查。SHA256SUMS.txt 校验本轮 C 交付文件。
