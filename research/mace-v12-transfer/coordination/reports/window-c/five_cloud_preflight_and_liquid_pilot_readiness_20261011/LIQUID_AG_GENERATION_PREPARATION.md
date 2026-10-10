# 六个 1250 K、64 Ag pilot 的生成与 QA 准备建议

这是准备设计；未生成坐标，未执行 MD/AIMD/DFT，不修改已归档任务、参数或数据角色。

## 已核对的来源

已只读检查现有 C checkout 中以下公开报告：

- `five_cloud_dft_parallel_plan_20261010/MASTER_TASK_MATRIX.csv` 和 `PILOT_RELEASE_GATES.md`，交接基准 `bb8ad7e4b725b634de2e2198d141df2bfe9d78b9`。
- `unified_potential_dft_data_plan_20261010/METHOD_AND_ROLE_RULES.md`。
- `unified_potential_independent_review_20261010/liquid_ag_pilot_plan.csv` 和 `EXTERNAL_REFERENCES.md`。

未访问或枚举封存结构及标签；未读取 A/B 活跃目录。

## 六项身份保持不变

| Task ID | owner | RNG seed | trajectory parent ID | 密度倍数 |
|---|---|---:|---|---:|
| C-D-AGL-1250K-rho097-seed01 | D | 125009701 | PF-AGL-1250K-rho097-seed01 | 0.97 |
| C-E-AGL-1250K-rho097-seed02 | E | 125009702 | PF-AGL-1250K-rho097-seed02 | 0.97 |
| C-E-AGL-1250K-rho100-seed01 | E | 125010001 | PF-AGL-1250K-rho100-seed01 | 1.00 |
| C-F-AGL-1250K-rho100-seed02 | F | 125010002 | PF-AGL-1250K-rho100-seed02 | 1.00 |
| C-F-AGL-1250K-rho103-seed01 | F | 125010301 | PF-AGL-1250K-rho103-seed01 | 1.03 |
| C-D-AGL-1250K-rho103-seed02 | D | 125010302 | PF-AGL-1250K-rho103-seed02 | 1.03 |

全部仍 `STRUCTURE_NOT_READY`；六项属于原 30 项子集，不是新增六项。

## 前置来源与盒子门

1. 必须先审定具有可复核方程/表号、适用温区、单位及不确定度的 `rho_ref(1250 K)`。旧稿的约 9.3 g/cm³ 是线索，不能替代经审定的参考密度。1250 K 接近银的实验熔点，且采样器自己的熔点可能偏移，不能仅用目标温度证明液态。
2. 为每项固定 `rho_i=f_i*rho_ref`、`N=64`、组成 Ag64、立方周期盒；若 rho 单位 g/cm³，`V[A³]=N*M_Ag/(N_A*rho_i)*1e24`，`L=V^(1/3)`，`PBC=(True,True,True)`。记录原子量来源、计算精度与实际写出的 cell。密度不确定度传播 `delta L/L≈(delta rho/rho)/3`。这里没有代入未审定密度或生成坐标。
3. 若输入采用非立方盒须另审体积、最短盒长、有限尺寸和 k 网格；不得静默替换。实际 RDF 截止必须小于最短盒长一半。
4. 采样器需要源网址/版本/许可、权重或势参数完整 SHA、元素支持和实际 MD 代码版本。建议优先原始 MACE-MP-0b3-medium；Ag-only EAM 可作独立来源候选，但需核验具体文件的来源和许可，并注意模型可能不支持液体。V18 保持冻结，不能因已适合界面就默认适合液体。
5. MACE/EAM 只产生候选坐标与采样诊断；采样势能、力、扩散和压力不得写入 `REF_*`、不得作为 DFT 真值或凝固验证结论。输出到独立 `sampler_*` 字段/文件。

## 真实生成：拟议且需单独批准

为避免 1250 K 的过热晶体，每个 seed 先走独立的熔化—冷却来源路径，而不是只给室温 fcc 晶体加速度。

- 用户应审定每条轨迹的初始来源、盒子密度、加热温度及控温方案、采样器、时间步、总步数/墙钟上限。可供审定的 ASE NVT Langevin 初稿：1 fs 步长，固定 cell 与密度，独立速度/随机力 RNG；高温熔化候选段约 1600 K / 10 ps、冷却到 1250 K / 10 ps、1250 K 平衡约 20 ps、取样约 20 ps。六条初稿共 60 ps × 6 = 360000 步。**这些是预算起点，不是已经批准或证明能熔化/平衡的参数。**
- 如高温段未出现液体结构特征则停止对应来源路线，提出另一个温度/路径请求；不得私自升温或增加步数。初稿单轨迹最高可请求 100 ps，六条最高 600000 步；达到批准的步数或墙钟上限仍未通过 QA 就停止，不将未平衡结构包装成液体。
- 不从一个已选液体快照复制、缩放后换 seed 来声称六个独立父源；每条使用独立初始化/熔化/热化历史。若共同使用相同 fcc 原型、共同熔化快照或其他源几何，记录共享 `source_root_family`。任务表的六个 PF 标识只能表示轨迹子族，不能抹去共同父源；开发/测试隔离以来源根族和 lineage 为准。这不改变旧 ID，也不把当前拟议训练候选改成测试。
- 固定每个 seed 的 RNG 算法和版本；将初速度、控温随机流、帧挑选种子以明确的子流派生规则记录。更换 seed 会产生新的来源记录；不得复用前一条轨迹作为另一条的独立起点。
- 不从较大的非周期局部切片直接裁出 64 原子再附加 PBC；这会制造边界/密度和邻居伪影。可建议未来独立授权 128/256 原子采样用于有限尺寸相态检查，但目前 64 原子任务应直接生成与保存完整周期盒；更大盒检查不是这六项自带权限。

## QA 证据和拟注册门槛

接受门槛需要在打开 DFT 标签前冻结；不能根据 DFT 误差回头选择容易评分的快照。以下是具体的准备规则，不宣称已经验证适用于本项目。

| 项目 | 必须保存和计算的证据 | 拟议判断/停止规则 |
|---|---|---|
| 结构身份 | Ag64、稳定 atom ID、坐标、cell、PBC、frame index、最小周期距离、密度、源 lineage、sampler hash | 组成/PBC/顺序与任务不符、NaN/Inf、重叠或超出审定采样器安全范围则停止；不能自行移原子修复后继续 |
| 稳定取样区 | T、sampler energy、压力/应力、体积、质心速度的时间序列和块平均；时间自相关与有效样本量 | 密度固定不能作为液态证据。目标密度各有不同压力，不能强制所有密度零压。取样前后的分块温度均值须与1250K在预注册容差/置信区间相容，T/energy/pressure无持续漂移；无稳定区则不选帧 |
| RDF/CN | g(r)计数、归一化、误差、bin宽、有效采样窗口；限制r<Lmin/2；第一峰及第一极小值；`CN=4*pi*n*integral(g(r)*r²dr)`与逐原子CN分布 | 需要宽化第一峰与短程有序、可识别第一极小值、没有持续晶体壳层模式；同密度两轨迹RDF应在块误差内一致。记录不同rho的配位变化，不用一个未校准CN固定值强行判液 |
| 局部晶序 | 局部q6、邻域平均q6、相干q6键/晶体连接团、必要的CNA/PTM；全局Q6只作辅证 | 用同箱/同密度独立固体和高温液体采样参考校准分布及相干键规则；没有通用q6数字门槛。持续的大型相干晶团、fcc/hcp特征保留或增长则拒绝“液体”称谓；全局Q6低不能排除晶粒抵消 |
| S(k) | 周期盒允许的波矢网格、方向/壳平均和采样误差 | 观察是否持续存在fcc Bragg样特征，结合q6/RDF/MSD；64atom的少量k模式不能证明长程无序或热力学熔点 |
| MSD | unwrapped坐标，扣除质心漂移，多时间起点，3D MSD(t)、误差、拟合窗口及速度自相关 | 排除弹道/笼困早期区后应出现可复现的长时增长，不是晶体振动平台。只用末端一个MSD值不能判液；以样本窗口内两个预注册滞后区间检查线性/斜率稳定，窗口不足则延长需审批或记QA未通过。64atom扩散值不能直接作宏观物性证据 |
| 独立采样 | 六条trajectory source SHA、初始化和速度记录、不同RNG子流；各rho两seed的统计结果、descriptor自相关 | 分别熔化/平衡；不同seed不等于不同结构来源父族。RDF/CN不相容或同一轨迹重复被当成两个seed则停止并查原因 |
| 快照多样性 | 帧间SOAP/邻居距离/CN/q6摘要，RDF差异及局部环境聚类；选择时间规则 | 一个稳定取样区选一帧，采样间隔至少超过相关时间；跨seed必要时max-min选代表帧，但选择准则冻结且不看DFT标签。六个快照只是最小覆盖探针，不能证明标签充分 |
| 去重 | 原始文件SHA、含组成/cell/PBC的精确几何哈希、周期平移/原子置换后的近重复检验；来源族账本 | SHA不同不说明几何不同；要排除包裹坐标、平移、重排和同轨迹近邻快照。descriptor仅做候选近重复检索，不能替代最终几何核验；不得接触封存坐标做去重 |

局部q6常用方法来源：Steinhardt、Nelson、Ronchetti, Phys. Rev. B 28, 784 (1983), DOI `10.1103/PhysRevB.28.784`；邻域平均序参量：Lechner、Dellago, J. Chem. Phys. 129, 114707 (2008), DOI `10.1063/1.2977970`。这里引用方法，不宣称已下载或在本轮计算。

## 持久归档与 SHA

每条轨迹在自己的 `DFT_WORK_ROOT/<task_id>/structure_generation/` 保存：

- `source_manifest.json`：任务/owner、source_root_family、trajectory parent、所有源 URL/许可与源 SHA、seed/RNG/substream、sampler代码/权重/势文件SHA、T_schedule、density_source、rho/cell、采样/选择规则、时间戳。
- 冻结 `generation_config.json`；初始结构、初速度、轨迹、unwrapped/image-counter 数据、thermostat/温度能量应力日志、step/time索引；MD restart 位于独立路径，不能和DFT checkpoint混用。
- `qa_summary.json` 和 T、RDF、CN、q6、晶团、MSD/S(k)、误差/相关时间的 CSV/图；记录PASS/FAIL/UNVERIFIED原因，QA脚本及版本SHA。
- 选定 `candidate.extxyz` 及 frame index/时间，纯几何字段为主；sampler数值独立保存，不能命名为 DFT label；`candidate_manifest.json`记录atom ID映射、cell/PBC、几何来源与所有精确SHA。
- write/read roundtrip后验证组成、顺序、原子ID、cell、PBC、坐标；对每个原始文件及最终输入生成 `SHA256SUMS.txt`。冻结候选后不得改变坐标、cell、原子顺序或许可来源；变化需重新审查并生成新版本SHA。

建议只把轻量输入、manifest、QA摘要和SHA提交各自分支；完整轨迹/restart放持久工作区及独立备份。GitHub只是协调账本，不能替代大文件和checkpoint的实际持久存储。

## 采样预算与停止条件

- 六条初稿360000力调用，最高拟议600000力调用。没有合法可复用的本项目液态64atom CPU采样实测，故目前墙钟为 `UNVERIFIED`。未来用户生成授权应先含首条100–1000步短采样测量（这也是MD，需要批准），记录seconds/step、CPU限额和峰值RSS，再由C重估。设实测单步s秒，墙钟小时 = 总步数*s/3600；不得引用DFT单点3–8h作为MD预算。
- 每10fs保存轨迹一帧时，初稿约36000帧，最大60000帧；Ag64仅float64 positions/velocities每帧3072bytes，合计约111MB/184MB，尚未包括cell、原子ID、时间、采样结果、索引、文件开销和restart。预留1–2GB总采样档案只是工程初稿，不是实际checkpoint测量；完整预算须结合输出格式、保存频率、备份策略和实际free bytes复审，避免逐步保存全波函数或不设上限日志。
- 重叠/非有限值、突然失温/压、无法熔化、反复结晶、时间步不稳定、来源/许可/hash不一致、数据角色冲突、持久盘/内存不足、批准步数/墙钟/盘使用超限时，停止对应生成任务并保留日志；不自动换势、升温、缩dt、延长或重启。

## 结构生成与 DFT 分開授权

### G-STRUCTURE 请求

明确列出六个task、采样器源/许可/hash、已核实rho和六个cell、每条完整温度路径与MD参数、总力调用及实测资源预算、QA阈值、输出上限与停止条件。用户批准后只可运行指定生成和QA，不可执行DFT。这里没有获得这一授权。

### G-NUMERICAL 请求

生成QA全部通过后，提供真实candidate SHA、父族、输入原子顺序/cell/PBC、DFT输入与PAW SHA；用户另批代表同一快照Gamma与2x2x2的k网格比较。如果不通过，不自动增密；申请3x3x3/更密的指定补充数值任务、预算和停止规则。各密度相差3%的网格有效性需明确覆盖极端盒；可先做一个代表，但不能无证据一律宣告所有密度通过。数值对照保持numerical-control用途，不能多算成训练构型或新父族。

### G-SIX-LABEL 请求

C审查固定几何数值对照后，单独列六个正式单点的最终k网格、hash、同一能量—力口径、主机资源和预算。GPAW26.7/PAW1.2.1/PBE/PW500eV/sigma0.10eV/convergence1e-5保持审定基线，`native energy`与force-consistent`free_energy`同时保存，并保存全部3D原子力。预计六个标签18–48h四rank串行墙钟，加已有计划mesh控制约3–16h；合计21–64h、84–256rank-hours仍是旧粗估，不是新实测。首例超过8h不排入后续，单例12h或160SCF上限仅为待批准门槛；任何SCF未收敛、NaN/Inf、hash/参数不符、资源越界、checkpoint身份不一致均停止对应任务、不自动改参数或重跑。批准pilot不放行剩余24项。

## 当前判定

`STRUCTURE_GENERATION_READY=NO`：待审定密度/采样器/许可、具体路径及授权、资源计量；本轮没有生成结构。

`DFT_PILOT_READY=NO`：没有六个真实QA合格且冻结SHA的输入；k点对照、持久盘和实际资源预算仍需放行。即使后来这些门槛满足，也必须等待相应用户计算授权。

64原子、六个DFT单点只检验局部高温液态环境的可获得性和数值成本；不会证明真实熔点、扩散、凝固速率、相变能力或熔坑物理可信度。
