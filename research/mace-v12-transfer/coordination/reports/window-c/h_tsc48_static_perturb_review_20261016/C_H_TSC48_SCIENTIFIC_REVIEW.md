# C — H 48原子 Ti₃SiC₂ 静态小位移候选审核

日期：2026-10-16（Asia/Shanghai）。审核基准 H commit：`e923008453d71077421b651cf19c6705420d2356`。

```
C_SCIENTIFIC_REVIEW = NEEDS_REVISION
DATA_ROLE = HOLD
KGRID_CANDIDATE = REVISE
RESOURCE_READY = NO
DFT_EXECUTION_READY = NO
execution_authorized = false
```

这是正式审核记录，不是执行授权。结构几何和补充计算的科学价值通过审查，但 H 报告对 k 网格中心方式的解释必须更正。37 项原始证据读取记录及 SHA 见 SOURCE_SHA256.csv，数值结果见 metrics.json；原报告、未改动的候选和母体等保存于 evidence/。没有运行 DFT、MD、AIMD、TTM、松弛或训练。

## 1. 结构身份及来源

候选文件 `TSC48_static_perturb_std0p025_seed2026101603_UNLABELED.extxyz` 原始 SHA-256：
`da3820c477023e0bb423f87550ed82ae0c2a571f53d294c055f78d7d052d5e8b`，2572 bytes。

从 H commit 读取的字节与报告所列生成提交 `cc65a3e6b75be28c06fda8608b4064662d714117` 相同。源说明列出的 `4ee3daab860a846fec4ef0c0ca4b031a36cfc87e` 母体也与 H 提交内母体相同；母体 SHA 为 `157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4`。

独立解析确认48原子、Ti24Si8C16、parent_row_id=1–48、元素顺序与母体相同。晶胞行向量（Å）：

```
 6.1150000000   0.0000000000    0
-3.0575000000   5.2957453441    0
 0              0             17.6235000000
PBC = T T T
```

周期最小距离为2.001610818 Å（输入行号6与24，Ti–C）。相对母体的单原子位移向量 RMS=0.041608321 Å，最大值0.072299923 Å；全体平均位移近零。Ti/Si/C位移 RMS 分别为0.0412071/0.0443148/0.0408022 Å。数值与每轴标准差0.025 Å的静态小扰动描述相容，但不重新声称复现了原始随机生成程序。输入没有DFT能量、力或应力标签，也不是高温轨迹。

## 2. 已有五种几何、价值及父族隔离

逐项检查已公开 baseline k10×10×4、Ti ID3 z±0.02 Å、C ID9 z±0.02 Å 的原始结构/能量力输出、summary和verification。五项均收敛、verification=PASS，输出文件哈希符合各自verification。baseline的多种k点计算只计为一个独特几何，不计作多份独立结构。

候选与五项均不是1e-6 Å容差下的逐行周期几何重复；比较距离见 metrics.json。它同时扰动Ti、Si、C及多个方向，可补充现有单Ti/单C z方向响应缺少的耦合局域受力、尤其Si和横向响应。**值得作为一项受控补充单点投入资源**，但一帧不能建立响应曲线、应力/弹性覆盖、独立泛化或高温稳定性。仍缺应变响应、更多独立父结构、表面/界面和真正热采样标签；它不能替代这些数据。

父族必须沿用 `TSC_BULK_PARENT_OR_CONTROL_FAMILY`（历史别名 `bulk_parent_48` / `V19_Ti3SiC2_48ATOM_BULK_CONTROL_FAMILY`）。新seed、全原子微扰和重新编号都不产生独立父族。当前E2原始角色冻结为diagnostic-only；本轮不修改旧标签，故DATA_ROLE=HOLD。

若未来申请将这一父族作为TRAIN，须另行批准全族用途登记，并排除同族结构作为独立DEV/TEST。现有numerical-control档案可保留原身份，通过独立用途映射记录获准训练用途；不得默默重写档案角色。DEV和TEST须来自来源独立、经结构等价/近重复审查的父族，TEST由独立evaluator管理。本轮不读取或枚举封存数据。

## 3. 网格：密度有依据，Γ居中解释不成立

母体8×8×4与10×10×4比较的两端原始输出、日志、summary、progress及verification哈希均已复核。native能量差−0.03512095 meV/atom、free能量差−0.01668959 meV/atom、力向量RMS差0.004617983 eV/Å，符合该母体的2 meV/atom与0.01 eV/Å门槛。它支持10×10×4作为小扰动结构的**初始网格密度**，不证明扰动结构自身k点收敛，也没有提供应力收敛结论。

关键差异：原始日志为 `MonkhorstPack(size=[10,10,4])`、MP shift=[0,0,0]；其实际坐标含x/y=±0.05、±0.15…以及z=±0.125、±0.375。这是偶数标准MP网格，不含Γ点。**MP shift=0不等于显式Gamma-centered网格**。H报告第2和第6节将二者等同，不能直接用旧误差值为新的显式Γ居中网格背书。

最小修订：保留10×10×4大小，改为与母体完全相同的GPAW整数Monkhorst–Pack约定，无显式gamma override/额外offset，并让最终输入预检记录实际k点列表或其哈希。若确实要用显式Γ居中方式，必须明确其为不同采样，现有对照仅支持网格密度；相应数值可比性仍需另行规划和批准，不能宣称已验证。本轮KGRID_CANDIDATE=REVISE，不启动任何网格测试。

## 4. 方法、软件与输出合同

PBE/PAW、PW500 eV、非自旋极化、Fermi–Dirac 0.10 eV、energy/density/eigenstates各1e-5、Pulay beta=0.05/nmaxold=8/weight=100、maxiter=160，与母体归档相容。母体日志另明确PPCG niter=2，最终锁应补写，不依赖隐式默认值。固定离子、固定晶胞，禁止松弛。

Ti/Si/C PAW应分别绑定：
- Ti: bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f
- Si: 2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0
- C: fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983

读取H已归档版本、4-rank和动态链接证据确认GPAW26.7.0/ASE3.29.0/gpaw-data1.2.1/OpenMPI5.0.7以及四rank链路先前通过。采用4 ranks、OMP/BLAS各1线程与4核配额相符；这些是归档证据，不冒充C对H当前会话的实时探测或新MPI测试。无需重装软件；实际执行前核验版本、PAW文件字节、runtime和可用资源。

未来唯一单点应记录native energy与force-consistent free_energy各自口径、48×3完整力、stress（注明Voigt次序/单位/正负号）、SCF收敛及迭代数、输入/代码/配置/输出SHA。仅收敛并通过归档验收后才产生有效标签；失败输出保持FAILED/INCOMPLETE。

## 5. 资源与正式执行边界

H报告的快照为4核配额、32 GiB内存限额、约28.16 GiB内存余量、28,828,061,696 bytes磁盘余量。C未远程测量H当前峰值内存和磁盘。28.83 GB比30 GB建议值少约1.17 GB；30 GB不是GPAW固定最低要求，不能仅据这1.17 GB差额否定单点。当前也没有足以批准例外的scratch/波函数checkpoint峰值、累计归档量和保留余量预算。

母体同网格耗时约3.26 h只是历史参考。随机扰动可能改变对称性、有效k点数和SCF难度。H提出至少5 h配额仍需平台实际证据；shell unlimited不能证明外部墙钟。checkpoint持久性和恢复、超时/写盘失败行为均未确认，故RESOURCE_READY=NO。

资源放行所需最小补证：明确外部墙钟配额及有界runner；量化内存/磁盘/checkpoint预算；验证独立工作根和持久备份/恢复通道。可审查专用持久卷或经验证的远端归档方案，不把某个挂载路径当成唯一实现。普通Git单文件超过50 MiB应先落实合适通道，不强推大型checkpoint；不得删除历史科研数据满足门槛。

本轮只需要纠正网格中心方式并补齐执行包和资源证据，不要求重新生成候选或重算已有五点。科学审核修订通过、资源全部就绪后，H仍须获得用户针对唯一单点的另行明确授权。当前不签发执行授权。
