# WMH 2006 纯 Ag 采样器：原文与实际参数来源复核

`SOURCE_PASS=YES`，`SOURCE_DOCUMENTATION_PASS=YES`。文件身份见SOURCE_RELEASE.json，独立哈希／参数记录见SAMPLER_PARAMETER_VERIFICATION.json。只读取用户提供原件、解析网格与表函数；未生成原子构型，未运行能量／力／MD／DFT。

## 原文与发布说明

Williams、Mishin、Hamilton，Modelling Simul. Mater. Sci. Eng.14(2006)817–833，DOI10.1088/0965-0393/14/5/002。用户论文PDF SHA为2a1e8ec6d59e21275e671bdfd1b54b9ebcc04ab85439954a0f7aa6612a23bd57；Table1解析函数cutoff5.99501 Å，第821页预测Ag熔点1267±5 K（实验对照1235 K）。

NIST发布说明2009-02-10，用户PDF SHA为262d9102f5a77572792b9a1d72f6c55c01372ec798f3321285d76246446a8c53，确认纯Ag参数名Ag.eam.alloy、setfl格式，2009-02-04由C.A.Becker从Y.Mishin三PLT经cubic-spline转换。理想fcc a=4.09 Å总能−11.3999998989 eV／4原子，即−2.849999974725 eV/atom，是原文经典势软件核验值，C没有运行该原子能量检查，不是DFT标签。

## 实际文件与校验

用户提供下载记录、SHA清单、Ag.eam.alloy及三个作者原PLT。C重算四文件SHA／字节数，与两份用户记录均一致。下载来源是用户提供的NIST公开URL取得记录，C没有宣称从远端服务器再次取证。Ag.eam.alloy：780452 bytes，SHA958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447。

ASE3.29原样解析：纯Ag、Z47、nr=nrho=10000、所有表值有限。文件mass107.868 amu；密度立方盒采用ASE默认Ag质量107.8682 amu，采样也使用后者，不以header悄悄改质量。实际cutoff5.995011000293092 Å是论文值的转换精度；不缩短cutoff、不重表、不由PLT重新转码替换官方setfl。

## 径向表尾段

公布文件cutoff比实际最后网格r[-1]=5.994411499193063 Å多一dr=.0005995011000293093 Å。C显式审核这一已发表setfl约定，11点尾段纯表函数检查：max|density|≈1.979e-14，max|phi|≈2.066e-13 eV。独立SciPy重算吻合。这是数值近零，不是数学严格零。

允许原ASE样条沿该公布文件单dr尾段到原header cutoff，policy为C_APPROVED_PUBLISHED_SETFL_ONE_BIN_ZERO_TAIL；代码必须同时检查一dr范围、实际尾段值及冻结SHA。不是允许任意文件宽区间外推，也不剪短势范围。这些是表参数评估，没有原子总能／力／候选结构。

## 科学使用与再分发

论文©2006 IOP。发布说明没有声明参数MIT／GPL／CC0或其它OSI许可，不能杜撰，也不把ASE软件许可等同于参数许可。但NIST明确提供并说明科学使用方法，用户提供从该来源合法取得原件并要求用于本研究，本项目科学使用依据已审核。potential.license_verified=true仅表示C已审核此科学使用依据，不表示OSI或再分发许可确认。

redistribution_status=UNVERIFIED_DO_NOT_REDISTRIBUTE_PARAMETER_FILE。GitHub只存元数据／SHA／引用／代码，不分发文章全文、setfl或PLT本体；D/E/F本地合法持有同SHA参数。来源URL见SOURCE_RELEASE.json。无需用户再下载文献或参数。

## 1250 K判读、PBC与有限尺寸

该势文献熔点1267±5 K，因此1250 K为近熔点／轻微过冷液体候选，不能预称为此势零压平衡液体。论文大胞共存熔点不能当作64原子NVT小盒相边界；实验1250 K密度及±3%用于候选环境覆盖，不是势EOS验证。

从1600 K混合／失晶序后降至1250 K，在预先冻结的生产窗检查RDF、Q6／晶序簇、MSD、T/E/P稳定及独立性。若晶化，不挑晶化前短帧冒充液态通过；应报失败／未建立液态窗，需修改协议时按既定门槛报告。60 ps／64原子与两个seed都不保证相态充分或防止结晶。

rc大于L/2但小于完整L。ASE EAM的PrimitiveNeighborList用整数周期offset并保留同索引不同像，支持这类多周期像；不能改成最小像／去重索引。RDF只到半盒高；没有有限尺寸、熔点、扩散或凝固验证通过的结论。

六原task ID／seed／trajectory_parent保持，共同参数／随机周期pack／melt-cool协议注册为同source campaign；不同seed不能充独立测试父族。各机复用已有ASE，不需LAMMPS或新安装。源门通过后仍等待owner独立持久工作根、预算与暖机签署；C不在本机替代D/E/F生成六结构。
