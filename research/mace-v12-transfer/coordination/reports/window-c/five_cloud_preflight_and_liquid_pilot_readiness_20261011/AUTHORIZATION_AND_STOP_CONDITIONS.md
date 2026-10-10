# 放行授权、资源和停止条件

## 当前状态

仅离线准备已获授权。六条结构生成、DFT数值试算、六个最终标签均NOT_GRANTED。LIQUID_AG_PILOT_RELEASE.csv逐项HOLD；G/H不执行这六项。原30分配保持，余24项继续HOLD。

## G-STRUCTURE：单独生成授权

申请须指定六个task ID、已核验rho(1250K)/六个cell、采样代码/模型许可SHA、初始化/源父族、RNG子流、温度路径、timestep/thermostat、最大steps/重试、QA阈值、输出预算和持久工作根。来源、resources/runtime/持久性通过前不申请实际运行。生成许可不包括DFT。

待审初稿1fs、每轨迹60ps（高温熔化10ps、冷却10ps、1250K平衡20ps、生产20ps），六条360000步；拟议上限100ps每条共600000步。是否真正熔化/平衡需QA证据，固定数量不是保证。实际CPU seconds/step与wall time UNVERIFIED；运行授权可先含首条100–1000步warm-up测量，它仍是MD，须明示许可。

资源容量规划：4CPU host、16GiB可用RAM建议不变；峰值RSS未测。时间预算为总steps乘实测seconds/step。例若0.1/1.0秒每步，六条60ps分别10/100串行CPU墙钟小时，D/E/F均分约3.3/33.3h；例值不表示实际性能。每10fs保存positions/velocities的六条基底数据约111MB/184MB（60/100ps），总工程预留1–2GB加备份/元数据并复核实际格式。

停止：来源/hash/许可/rho不一致；重叠/NaN/Inf、失温压或持续漂移、无法熔化/重新结晶、MSD平台或窗口不足、独立性/重复问题、内存/盘/steps/wall budget越界。停止对应新任务并保存证据；不自动换势、升温、减dt、加时/改seed或重跑，不影响既有任务。

## G-NUMERICAL：单独DFT数值授权

真实候选QA/源父族/结构SHA冻结后，用户另批同一代表snapshot的Gamma与2×2×2比较；3×3×3或更密必须另批，并重新审盘/内存。各密度极端盒的有效性仍需审查，不以单例证明全部。初始1–2项mesh控制仍只是预算建议（约额外3–16h/4ranks），不能自动运行。

方法维持GPAW26.7、ASE3.29、PAW1.2.1、PBE/PW500eV、Fermi sigma0.10eV、SCF density/energy/eigenstates各1e-5的既有GPAW输入语义；maxiter160为拟议上限，须在计算授权中冻结。same-geometry mesh comparison待审目标deltaE≤2meV/atom、force RMS delta≤0.01eV/Å、max force delta≤0.03eV/Å；不通过即HOLD，不静默改科学设置。数值子运行保留数值控制角色。

输入冻结组成、atom order/IDs、cell、PBC、kmesh、PAW/方法SHA及checkpoint身份。保存native energy (force_consistent=False)、force-consistent free_energy及完整3Dforces；热标签接入统一free_energy+forces口径，不能混不同energy keys。

## G-SIX-LABEL：六个最终单点授权

C审查SCF/mesh/force/hash、RSS/checkpoint实测及资源预算后，用户单独批准六个输入SHA与最终k点。旧粗估每项4MPI ranks 3–8h、六项18–48串行wall h=72–192rankh；加数值控制共21–64h=84–256rankh。各D/E/F有2label任务、各6–16h，mesh控制和主机速度另计。不是实测承诺。

停止：SCF max160不收敛、NaN/Inf、身份/hash/checkpoint/参数差异、重复task、能量力口径不一、4CPU/memory限额不符、预测写盘越过风险边界、OOM风险或超批准wall budget。待审首例8h复核（不启动后续）、单例12h hard cap；必须与授权单一致。不自动重启/调参/清理。采用已批准安全flush方式停对应新任务，保留最近完整checkpoint与失败日志，禁止涉及A/B或他人进程。

## 磁盘例外和更密网格

30GB原建议不变。五台现有available均低于此数值；若请求有限Gamma/2mesh例外，须先有真实bands/IBZ/FFT与checkpoint测量/保留账本、明确临时写盘副本/风险余量、可靠持久卷和user专项批准。本轮不改变门槛。更密mesh可需要31–50GB甚至更多，即便达到30GB也不自动满足。

15Ag/10TSC旧numerical controls不转训练；历史评分父族不转独立测试；不读取或枚举sealed数据。静态标签合格不能等同熔化、凝固或熔坑适用性通过。
