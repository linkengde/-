# 结构采样与DFT写盘的独立放行

30GB是原建议，不是GPAW固定最低要求，保持建议但不作为小型几何采样禁止条件。当前owner28.8–29GB读数不能自动判MD不可行。也不能据此自动允许大型DFTcheckpoint。

## 六轨迹采样

每项Ag64，拟议120000步/.5fs=60ps；每100步保存coords/momenta与连续unwrapped positions，每20步热日志。两任务在各owner的独立sampling root；单轨迹档案cap512MiB（含restart bundle），开跑前可用量≥2GiB，RSS上限8GiB、warmup累计≤1h、continue累计≤5h、每项全程累计≤6h为待源审查后冻结的拟议上限，不保证相态平衡充分。需owner提供真实root路径、挂载/持久策略、剩余quota和相互隔离证据。不能用未验证overlay声明持久；若使用远程备份，须给明确上传/校验频率、恢复方案和失效风险供C审核。当前未创建远端root。

warmup1000步包括在hot20000中；先C通过来源/rootgate，再owner跑已批准六轨迹内的短段，回报seconds/step、RSS、输出字节、T/E/P稳定及restart身份；C逐task加入continuation_task_ids才继续。此短段并非独立新构型/新增DFT任务。原seed/PF不改，真实source campaign已登记，六条实际轨迹尚未生成。

## GPAW checkpoint

先真实候选通过QA和SHA冻结。D Gamma/2MP compare在同几何；保存actual bands/IBZ/PW/FFT/RSS/checkpoint大小，单active SCF/host。current+previous+atomic temporary与已有mesh文件分别登记；至少5GB风险reserve，并由C独立run config绑定实际preflight/input/方法hash审定capacity。原28.8GB不直接判不够，预算证据不够时仍HOLD；不清理旧checkpoint。

3mesh、更密、超12h/160SCF、改cutoff/sigma/SCF科学参数均需用户专项批准；Gamma/2不通过则暂停，不擅自增密或偷用Gamma。三密度有效性仍须核对，不能一个代表默认全部PASS。

## 当前资源资料需求

D/E/F分别提供自己的provider持久卷说明或经验证备份/恢复方案、sampling root和DFT_WORK_ROOT、time-stamped byte/quota读数、runtime SHA。不重复安装/全套软件预检，不增加云端，不改30任务总分配。E补memory.current，纠正system available超过32GiB限额的口径。

SOURCE_GATE已true；STORAGE／WARMUP／CONTINUATION_GATE在配置中目前false；普通owner自报true不会自动解除C签字门槛，发布config SHA版本后才可执行。
