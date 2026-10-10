# 五云端预检汇总与首批液态 Ag pilot 准备

## 结论

PREPARATION_DOCUMENT_READY=YES；HOST_READINESS_VERIFIED=NO；STRUCTURE_READY=NO；PILOT_RELEASE=HOLD。

本轮未执行MD/AIMD、DFT、训练、LAMMPS或TTM，未生成坐标。A/B的实际工作区、进程、input/checkpoint/参数/任务卡不变；C checkout仅纳入已有公开main历史，无A/B文件内容修改。

用户确认五台安装齐备和资源独立，本次五台owner实时反馈已全部收到。无D–H原始已归档报告，也无C远程连接，所以正式 [HOST_READINESS_MATRIX.csv](HOST_READINESS_MATRIX.csv) 区分OWNER_REPORTED、USER_REPORTED和UNVERIFIED，不将任何单机数据外推给另一台。实时反馈齐全不表示证据完备；五台底层资源独立性仍需provider证据。

## 主机汇总

| host | 实际CPU额度 | owner提供可用内存 | available disk bytes | runtime重点 |
| --- | --- | --- | --- | --- |
| D | 4，quota400000/100000 | 30,804,070,400 bytes，限额32GiB | 28,797,075,456 | onboarding/runtime.sh语法通过；cgpaw backend下4rank成功，裸Python探针失败 |
| E | 4，quota400000/100000 | 报34,777,948,160 bytes，大于32GiB限额，effective余量未确认 | 28,800,167,936 | .setup/e/runtime.sh语法/setup PASS，work-root BLOCKED |
| F | 4等效cgroup quota | 约28.64GiB cgroup余量，需精确bytes | 28,874,539,008 | .setup/onboarding/runtime.sh可加载；setup/library/historical MPI PASS；work-root unset |
| G | 4，quota400000/100000 | 约28.9GiB，限额32GiB，需精确bytes | 29,018,927,104 | runtime.sh缺失；已有activate.sh此前setup解析匹配 |
| H | 4等效，需原始quota回执 | 34,204,598,272 bytes，需口径/cgroup补证 | 28,896,153,600 | .setup/runtime.sh；setup/MPI linkage owner报通过，需SHA/完整回执 |

可见5vCPU/cpuset不等于5核计算额度。上述数值来自owner消息，非已归档主机报告或C实测。五台持久性/配额/独立性仍未验证；D/G/H工作盘为overlay或owner未核验持久性，不自动作为持久DFT卷。

F/G/H沿用已有软件，不重装。F的DFT核心栈owner已核验；缺LAMMPS不阻碍当前GPAW准备，本阶段无需安装。其ENV_READY=NO是owner更广范围状态，不意味着DFT核心栈需要重装。E所报available大于32GiB限额，需要区分系统MemAvailable与cgroup剩余，本报告不据此认定有效可用额度。D的cgpaw backend实测要求为其他同build主机的核验线索；G缺runtime属于配置准备缺项，不能据此判需要安装。详见 [WORK_ROOT_AND_RUNTIME.md](WORK_ROOT_AND_RUNTIME.md) 与 [READ_ONLY_PREFLIGHT.md](READ_ONLY_PREFLIGHT.md)。

## 磁盘审查

30GB建议不变，本轮不批准降低或例外。D/E/F/G/H分别约28.797、28.800、28.875、29.019、28.896十进制GB，距30GB差约1.203、1.200、1.125、0.981、1.104GB；精确字节与GiB在矩阵中。还需卷配额、记录时间和累计保留账本。

公开48atom TSC历史checkpoint为7,361,378,540 bytes；它说明compact结果很小不等于restart很小，不能直接套算Ag64。64Ag500eV的bands/IBZ/PW/FFT尚未实测。保守三份checkpoint（current+previous+atomic temp）、已有mesh文件、日志与风险余量的工程示例给有限2mesh约16.8–23.4GB，3mesh约31.0–50.4GB；假设rho=9.3、bands662和overhead均未验证，不能据此声称当前磁盘够用。9.3仅用资源尺度示例，不用来生成坐标。详见 [DISK_AND_RESOURCE_BUDGET.md](DISK_AND_RESOURCE_BUDGET.md)。

## 独立目录与Git

[WORK_ROOT_PLAN.csv](WORK_ROOT_PLAN.csv) 提议各自仓库外持久根 /workspace/dft-work/cloud-x/ag-ti-si-c-b1；不使用/tmp、共享checkpoint或A/B目录。仅命名隔离仍需挂载/平台分配证据。没有远端创建/配置操作。

handoff bb8ad7e4b725b634de2e2198d141df2bfe9d78b9 未进入所审查main 684e9a887547404d959e2a7b679c9596ea134a8a，两者分叉于bcb51e5。C新专属分支保留旧handoff并普通merge main，无冲突，不push main。[BRANCH_AND_REGISTRATION.md](BRANCH_AND_REGISTRATION.md) 给安全合并/登记提案。D消息中的f954仓库SHA是历史快照定位，不是资源报告提交，不能当作D报告已归档。

## 液态 Ag 准备与授权

1250K接近常用熔点1234.93K；NIST原页面未取得，rho_ref(1250K)仍UNVERIFIED，不填最终六个cell/坐标。原始MACE-MP-0b3-medium的官方模型目录和MIT许可已获取，现有权重SHA一致，推荐作为几何采样候选；液态能力尚需QA。Ag EAM具体来源可核验，但文件许可适用范围/液态适用性仍待审，暂不选。详见 [SAMPLER_AND_DENSITY_SOURCES.md](SAMPLER_AND_DENSITY_SOURCES.md)。

[LIQUID_AG_PILOT_RELEASE.csv](LIQUID_AG_PILOT_RELEASE.csv) 保留六个旧task、seed和D/E/F各2项；是原30任务子集。每项STRUCTURE_NOT_READY、输入SHA未计算、生成与DFT授权均未授予。计划PF只是trajectory ID；另加source_root_family防止共同初始化后代冒充独立父族。原数据角色完整保留，不将候选改为训练/测试。

[LIQUID_AG_GENERATION_PREPARATION.md](LIQUID_AG_GENERATION_PREPARATION.md) 给真实熔化/冷却/平衡设计、RDF/CN、局部q6/晶团、S(k)、unwrapped去COM MSD、独立采样和精确/近重复检查。60ps/seed、六条360000步只是待审初稿，不保证足够。温度/低全局Q6或随机扰动不能单独证明液态；64atom也不能证明熔点/凝固物性。采样势energy/force与DFT标签分开。

## 推荐下一步与停止

1. 将五台owner反馈归档为原始报告，补时间戳/完整回执/配额/持久性、F/G/H fresh-shell runtime脚本SHA与代码审读，并核定E的effective memory；补rho1250文献表/公式与不确定度。
2. 通过后单独批准六条真实结构生成和QA，包括温度路径/预算；当前不生成。
3. 冻结QA合格六输入及SHA后，单独批准固定几何kmesh数值试算；C审核后用户再批准最终六个单点。
4. 余24项保持HOLD。若拟请求磁盘例外，先补实测/保留策略证据和风险评估，再由用户批准，不自行改30GB建议。

详见 [AUTHORIZATION_AND_STOP_CONDITIONS.md](AUTHORIZATION_AND_STOP_CONDITIONS.md)。当前最小建议是补来源/持久资源证据，科学计算放行仍NO。

更正旧文书术语：TSC10×10×4是48atom reference cell，不是primitive cell；只在本目录更正，不覆盖旧报告。15Ag/10TSC数值控制角色不变，不读取/枚举sealed标签或结构。交付SHA256SUMS覆盖全部新文件，本轮完成后停止。

发布前main更新到284285cec6a7a5691ffc01fd9616d7760dc4dac4，只有B公开归档/状态变化，C已普通merge同步且无冲突；细节见分支/来源记录。科学输入和资源审查不因Git同步改变，新增B标签未展开分析或变更角色。
