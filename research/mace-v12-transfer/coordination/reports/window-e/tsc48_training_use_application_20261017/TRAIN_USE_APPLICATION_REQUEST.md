# E → C：bulk_parent_48 六点有限探索性 TRAIN 用途申请

状态：**PENDING_C_EXPLICIT_USE_APPROVAL**。本文件提交审批请求；没有改变原始数据角色，没有创建 TRAIN 映射，也没有授权实际拟合或其他计算。

## 申请事项

请 C 审核是否允许将以下六个已完成的 DFT 标签单独登记为 `bulk_parent_48` 父族下的 **TRAIN 资格**，仅用于后续受约束、力重点的 Ti₃SiC₂ 骨架 TSCBOP 探索性修复准备：

- `TSC48_VOL_m02` — input SHA `1e8aeed059783a36e471873a3b3821369750677ae93ec5965d1d5760189ed683`
- `TSC48_VOL_p02` — input SHA `a86ee4560676e1a9679729cd6c9721b4b4a875d41cc80862d6287a056ddb0db1`
- `TSC48_BASAL_XY_SHEAR_m005` — input SHA `72be50f6e2fc3c7fcb610fadf520a7be3eadb8b14c2c7f27416da2fe7881ed7b`
- `TSC48_BASAL_XY_SHEAR_p005` — input SHA `cd0e3d09e413d25a54dbaa38868106934b71e2d95c3c879c8346dbd0cff8f3a1`
- `TSC48_SI_ID7_Z_m020` — input SHA `77ab3ae7dc59a41adea035003543fe19ae9c8deb8dd2e1413807ed25164e91f4`
- `TSC48_SI_ID7_Z_p020` — input SHA `8222f1538719c859765f1ab48db483e4599d4f64ac2835f67dc5709ee860002d`

六项属于同一个 48 原子母体结构族，不能拆分成相互独立的训练/验证/测试家族。此申请只含上述六点；不包括旧五点、H 候选或其他后代。

## 拟议用途边界

- 初始目标只使用逐原子力（eV/Å）；能量与应力保留为诊断，不进入初始目标函数。
- 不启用 stress loss，不据此标定弹性常数。
- 不把同族数据当作 DEV/TEST，也不声称独立验证、泛化成功或生产可用。
- C 若批准的是 TRAIN 用途映射，也不代表参数优化、DFT、MD、TTM 或 MACE 训练获批；任何实际拟合仍需单独的用户授权。

## 已知限制

- 实际网格为标准 Monkhorst-Pack 6×6×2。母体网格对照只覆盖未扰动母体，未证明六个扰动结构或应力收敛。
- PAW 身份证据为批次级哈希记录与一致性核验，不是每项作业/每个 MPI rank 的历史现场哈希。
- 五项 launcher 退出码仍未知；C 接受其原始 SCF 收敛与完整标签作为有限探索用途的来源证据。
- native energy 与 force-consistent free energy 为不同口径；本申请不将二者用于初始拟合目标。

## 审批后状态处理

在 C 明确批准前，`TRAIN=[]`、`DEV=[]`、`TEST=[]`，原 `DIAGNOSTIC_ONLY` 归档身份不变。获批后也只允许创建逐项、可追溯的用途映射；不覆盖原始 manifest。若 C 拒绝或要求更窄范围，保持 HOLD 并按意见修订。

C 的来源复核提交：`{c}`；E 补证提交：`fa699e1340e8777a7095bb8c639a7247c28d1956`；标签资格审核提交：`b3eef435e372b833053b369b4d401165be28f745`。
