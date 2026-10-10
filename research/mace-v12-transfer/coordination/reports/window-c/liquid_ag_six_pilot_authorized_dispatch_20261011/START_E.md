# 云端 E：本轮六液态 Ag pilot 启动交接

仓库linkengde/-。原交接固定bb8ad7e4b725b634de2e2198d141df2bfe9d78b9；C本轮branch cloud-c-liquid-ag-six-pilot-execution-20261011。先正常fetch并只读检查本包SHA；不覆盖活动input/checkpoint，不强推main。原HANDOFF中其余4项HOLD，不计算他人的task。

用户已批准以下两项的真实生成和合格后的DFTpilot；不需要再次请求整轮运行授权。C的科学/resource gates仍必须通过，本包目前SOURCE_RELEASE=YES，但STORAGE／WARMUP_RELEASE=NO。不能把“读到启动文件”当成质量审核通过。

| task ID | planned seed | target rho g/cm³ | cubic edge Å |
| --- | --- | --- | --- |
| C-E-AGL-1250K-rho097-seed02 | 125009702 | 9.032733418760 | 10.826824658074 |
| C-E-AGL-1250K-rho100-seed01 | 125010001 | 9.312096308000 | 10.717455315296 |

每项Ag64，PBC true/true/true。密度及论文液态温区已由上传原文独立核验；1250K采用邻近表值保守工作误差0.26%，具体采样势文件和SHA已核验，当前没有六个坐标或可运行DFT输入。parent保持原trajectory身份，同时新增真实source_root lineage；同源seed不变成独立测试。

## 当前立即执行的范围：只读补证

独立核对D返回的势文件身份/SHA/许可/实际cutoff与ASE周期图；补memory.current和effective headroom，所报available大于32GiB限额不可作cgroup余量。
各自执行preflight_owner.py只读采集新资源/环境，提交本机报告、时间戳与文件SHA到自身branch。effective4CPU、可用RAM、work-root持久性/隔离、实际remaining bytes和checkpoint retention ledger必须明确。helper输出是本机自报，不是C远程验证。不要安装LAMMPS或盲目运行安装脚本/H2 smoke。

论文／发布说明及实际 `Ag.eam.alloy` 均已核验；不能把1986 Ag_u3.eam改名替代，不能擅自换MACE/EAM。把原文件SHA、取得来源、发布说明引用交C核验；项目科学使用依据和再分发权限分开记录，未知再分发许可时不上传参数本体。密度原文已齐。若合法源访问也403，记录失败并报告，不绕过。文件证据轻量提交自己的branch，full轨迹/checkpoint将来放独立持久根。

## C通过source/parameters gate后

只用C冻结的具体sampler文件SHA、cutoff、密度、初始化和温度路径、time step/thermostat/时长/seed。先真实消除晶序、1250K定密度平衡再抽帧，不将静态高温或同结构换seed伪称液体。不随意缩短、升温或加时。每项完整轨迹；RDF/CN、局部q6/CNA、MSD、T/E/P时间序列、周期距离与重复度通过后送C。未经C确认不得DFT。

`run_sampling.py` 已提供；当前模板仅提案：.5 fs、1600 K混合10 ps／冷却10 ps／1250 K平衡20 ps／生产20 ps，暖机1000步包含在总120000步内。实际参数与QA门槛尚未签署，不凭模板启动MD。该势文献熔点1267±5 K，1250 K属于轻微过冷候选，生产窗晶化必须报告失败，不挑晶化前帧替代。64Ag立方盒不能直接由整数fcc conventional重复得到64原子；初始化方法需独立说明，禁止裁剪大体系后补PBC制造边界伪影。

## C通过structure gate后

等待D的同几何mesh对照及C放行最终k点，不先跑自己的Gamma单点。
GPAW26.7/ASE3.29/PAW1.2.1，PBE PW500eV Fermi sigma0.10eV。SCF energy/density/eigenstates各1e-5沿用精确GPAW口径；max160和wall12h为用户本轮批准上限；首例>8h暂停后续并重估。同结构energy delta≤2meV/atom、force vector RMS delta≤.01eV/Å、max atom delta≤.03eV/Å；能量键明确。

保存native energy、force-consistent free_energy、完整3D forces、atom order/cell/PBC、方法/PAW/结构SHA、SCF/kmesh、runtime/RSS/disk/time/checkpoint状态。每task每mesh独立run目录；all-rank collective temp checkpoint再原子替换，不能直接套用B会自动改任务卡的runner。30GB建议保持，任何例外要实际预算/风险证据；不清理旧checkpoint来凑空间。

来源/相态/hash/SCF/finiteforce/资源/预算失败即隔离对应任务、保留证据，不自动改科学参数或重跑。结果按STRUCTURE_PASS、SCF_PASS、KPOINT_PASS、PROVENANCE_PASS、PILOT_ACCEPTED逐项回报。完成自己的2项即停止；G/H待命，其余24、训练/生产MD/TTM不启动；A/B不接管不改动。

## 程序使用顺序（只使用C发布的配置SHA）

以下为执行接口说明；当前模板放行门关闭。预检是只读，实际文件与运行目录由本机核实。

```bash
python PACKAGE/preflight_owner.py --owner E --work-root VERIFIED_DFT_ROOT --runtime VERIFIED_RUNTIME > NEW_E_REPORT.json
python PACKAGE/inspect_eam.py --potential LEGALLY_OBTAINED_Ag.eam.alloy --format alloy --output NEW_INSPECTION.json
python PACKAGE/run_sampling.py --owner E --task-id OWN_TASK_ID --config C_FROZEN_SAMPLING.json --config-sha256 C_PUBLISHED_SHA --stage warmup
```

C审查暖机后签署新的continuation配置SHA；同一任务使用`--stage continue`。每台两轨迹串行，单进程运行，不用MPI。输出根和参数文件路径必须来自C配置，禁止使用他机或A/B目录。

```bash
python PACKAGE/qa_sampling.py --owner E --task-id OWN_TASK_ID --task-root REAL_TASK_ROOT --config C_FROZEN_SAMPLING.json --out NEW_QA_DIRECTORY
python PACKAGE/run_dft_pilot.py --config C_APPROVED_DFT.json --config-sha256 C_PUBLISHED_SHA --owner E --launch --execute
```

DFT runner自行建立本机4-rank新作业组和12 h墙钟监督；不要再套外部MPI／setsid／timeout启动。关闭门检查在GPAW导入和SCF前退出。实际DFT config绑定C签署source/QA、inputSHA、机器fingerprint、runtime/PAW与工作根；模板不能自行改true。

C只汇总自己的新增报告。本人分支保存轻量来源／QA／结果清单；checkpoint、参数表和完整轨迹留独立持久根，未知参数再分发权不上传本体。

## 已冻结采样器

`Ag.eam.alloy` SHA-256：`958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`，780452 bytes，纯Ag setfl，文件cutoff=5.995011000293092 Å。来源：

https://www.ctcms.nist.gov/potentials/Download/2006--Williams-P-L-Mishin-Y-Hamilton-J-C--Ag/2/Ag.eam.alloy

本地下载必须校验此SHA。公开科学使用依据经C审核，再分发权未确认，因此不要将参数本体push到GitHub。不能以不同Ag势替换或按PLT重新转换。文件尾段policy和1250K轻微过冷判读见SOURCE_RELEASE.json与SAMPLER_SOURCE_REVIEW.md。当前只需补本机独立持久根、剩余容量及runtime回执，随后C发布暖机配置。
