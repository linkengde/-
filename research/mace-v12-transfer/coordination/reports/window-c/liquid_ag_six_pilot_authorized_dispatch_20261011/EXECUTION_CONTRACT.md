# 执行与接收记录的契约

## 来源和数据角色

SOURCE_RELEASE.json只签署资料／参数身份与本项目科学使用依据。实际文件纯Ag setfl，原样保留SHA与cutoff；不改权重、不重转PLT，不以经典势E/F当DFT真值。六项共享同一source campaign父族，保留各trajectory_parent／seed；生成方案不是已存在坐标，seed不同也不是独立测试来源。

15个Ag与10个Ti3SiC2旧数值控制及B四位移点角色保持；本六项先作为公开候选数据，不自动加入训练。不得读取封存结构／标签，不能用已评分父族充独立测试。

## 采样配置：提案到签署

C_RELEASE_CONFIG.json当前source=true，storage／warmup／continuation=false。实际采样路径及各机root尚未签署。后续C发布按owner绑定路径与证据的配置及SHA，不允许本机把false自行改true。

协议提案为单进程ASE3.29 EAM、64 Ag、固定周期立方盒，默认ASE Ag质量107.8682 amu。独立seed的PCG64子流分别用于周期随机pack、速度与thermostat；初始最小周期距离2.2 Å，不裁晶体、不挪盒边避免拥挤。Langevin .5 fs、friction .01/fs，1600 K 20000步（含1000暖机），降温20000步，1250 K平衡40000步，预定生产40000步。120000步／60 ps是预算，不保证足够。

暖机累计限1 h、续算累计限5 h、每项全程累计限6 h、采样RSS上限8 GiB、可用空间保留2 GiB；单轨迹档案含重启副本上限512 MiB。各owner两任务串行，每步资源门与显式停止；续算不得通过重新调用绕过累计墙钟限额。精确NPZ动量、PCG64 RNG、raw unwrapped坐标、全文件fsync与原子状态发布供重启检查。所有owner输出路径不得用symlink逃出自己的root。

QA提案已登记于config，但thresholds_frozen_by_C=false。通过暖机和来源审核再冻结，不能事后调整门槛迁就轨迹。必须检查所有归档phase最小距离／T，预定生产窗RDF/CN、局域q6／coherent bonds与晶序簇、T/E/P四块、去COM多起点MSD。晚时log-log指数与三窗口斜率一致性用于排除弹道或平台运动；这些仍是筛查规则，不能替代物性／凝固验证。

AUTO_DIAGNOSTIC_PASS只产未签pure-geometry candidate；STRUCTURE_PASS始终false。C另审高温失晶序、是否持续晶化、轨迹来源与六项重复度，再签署冻结输入。RDF仅到半盒高；共同父族和64原子有限尺寸限制保持。

## C签署的DFT接收文件

直接把qa_report或raw source_manifest改true不符合此契约。C在独立新接收目录另存manifest，并由DFT配置绑定这些文件SHA：

- 来源接收记录：SOURCE_PASS=true、真实source_manifest／trajectory／potential／configSHA、task／owner／seed／parent／source_root及明确的数据角色；原来源文件不变。
- 结构接收记录：task_id、owner、STRUCTURE_PASS=true、PROVENANCE_PASS=true、input_sha256（候选文件原字节SHA）、geometry_sha256（实际extxyz读回几何）、localQA与六项重复审查SHA及C结论。
- 主机接收记录：原只读owner_preflight及其SHA、C审核持久／隔离证据。runner所需owner／work_root／runtime_path／runtime_sha256／machine_fingerprint_sha256保持与本机一致，C在新接收记录签persistence_verified／independent_resources_verified；原自报不覆盖。machine fingerprint只绑定本机身份，不证明物理资源隔离。

纯几何input必须Ag64，atom_ids=1..64、原子顺序固定、cell/PBC匹配；不带sampler结果calculator。geometry SHA按numbers／12位positions与cell／pbc／atom_ids的compact sorted JSON计算，使用序列化读回值。

## DFT单点与停止

DFT_RUN_TEMPLATE.json保存PBE／PAW／PW500 eV／Fermi .10 eV、energy/density/eigenstates各1e-5原API口径、max160等方法。energy阈值是eV/valence-electron，density是electrons/valence-electron，eigenstates是eV²/valence-electron；不得将1e-5称为总能精度。native energy与force-consistent free_energy分别保存，forces全原子eV/Å。

已授权的D固定结构Gamma及标准shifted MP2×2×2先对照。双energy差≤2 meV/atom、全原子向量force-RMS≤.01 eV/Å、最大逐原子delta≤.03 eV/Å是数值筛查目标。C另审其它密度适用性；一个快照不会自动放行所有密度。3×3×3或更密需专项批准。

runner自行启动localhost四rank（不读取远程hostfile），MPI/BLAS每rank单线程，不套外部MPI／timeout。watchdog只管理自己新session，8 h暂停后续，12 h结束本作业并保留最近完成checkpoint；160 SCF仍不收敛则FAIL。失败／初始化异常一律暂停本owner后续，禁止自动重启、改参数或清理旧checkpoint。

原30 GB建议保持；只在C按实际bands／IBZ／PW／FFT／wavefunction dtype及当前磁盘quota审核新pilot有限预算后签低于建议的例外。current＋previous＋临时写盘和其它已保存mesh分别计入，至少5 GB reserve。脚本仅管理本新pilot自己的两代checkpoint；不碰A/B及别人的文件。

C可引用合格同几何／同方法的2网格control作为该pilot的公开参考候选，保留来源／数值控制角色，不重复相同单点，也不改变训练资格。six-pilot仍是六个计划构型而非“控制另加六个”。

## GitHub与完成判据

C提交代码／来源元数据与审查；owner在独立branch提交compact来源／QA／结果／SHA。实际参数、全轨迹与.gpw留在独立持久根，当前未知参数再分发许可时不push本体。正常fetch/commit/push，出现冲突停止，不强推main。

源通过不等于暖机通过；暖机通过不等于液态通过；数值收敛通过不等于训练资格或凝固性能通过。没有真实运行记录时，各科学结果均NOT_RUN，不能用脚本测试替代六项DFT结果。
