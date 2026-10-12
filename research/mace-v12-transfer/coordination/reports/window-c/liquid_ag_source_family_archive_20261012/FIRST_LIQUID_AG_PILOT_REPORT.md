# 首条液态 Ag 轨迹试运行报告

日期：2026-10-12  
节点：D  
目标：`AgL_1250K_rho1p00_seed01` / `SF_AGL_1250_R100_A`（计划 ID）  
状态：`PREFLIGHT_BLOCKED`；未生成结构，未启动 MD。

## 结论

本次没有启动试运行。现有方案和登记表将采样器、精确密度、随机初始化细节及积分/热浴参数列为待冻结或待同步内容。用户要求这些项目不能核实时停止、不猜测参数；因此没有生成坐标、没有计算力/能量、没有执行 dynamics step，也没有产生 trajectory、checkpoint 或伪造的运行统计。

## 已完成的设计归档

现有三份设计文档和候选登记 CSV 已复制归档，并在 D 独立分支上传。归档包内 `SHA256SUMS.txt` 对设计、登记表、README 和本报告逐文件校验。设计归档提交：`a2bbc97400ffa0a94493fb26c333e32691080700`。该提交经独立临时 Git 仓库 fresh-fetch；逐个读取远端对象后，所有文件 SHA-256 与清单相符。

归档文件身份：

- `LIQUID_AG_NEXT_PHASE_PLAN.md`: `47b95ecd33dccdb2429638c0bb511dc900eca7f50347b1923e75fc03b9e506aa`
- `LIQUID_AG_SAMPLING_DESIGN.md`: `c6619ac69fd1f9fe5de388986cd9a8d044f2edcd0af5fd1ec6d7afcf85bc9080`
- `LIQUID_AG_SOURCE_FAMILY_DESIGN.md`: `920fd43cdfc3bb0ec4b5c8ad831cc7a9421018fed70752193f14229774d62a2c`
- `liquid_ag_pilot_plan.csv`: `563b36d49e033843ab868edd50523297cc6c9c97d86bd37b255e5fd367300309`

## 实时 preflight 证据

| 检查项 | 结果 | 证据 |
|---|---|---|
| Python / ASE | 通过 | Python 3.12.14；ASE 3.29.0，解释器 `/workspace/.venvs/gpaw-mpi/bin/python` |
| MPI | 已识别 | Open MPI 5.0.7；当前任务未运行 MPI 计算 |
| LAMMPS sampler | 缺失 | `lmp`、`lmp_mpi` 均未在 PATH 中找到；方案推荐 LAMMPS，但没有另一套已批准的可用 sampler 配置 |
| 势文件身份 | SHA 通过 | `/workspace/WMH2006_Ag-local-validation/Ag.eam.alloy`，780452 字节，SHA-256 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447` |
| 势文件来源/许可 | 未核实 | 当前可见计划要求在执行前核验来源与许可；未找到随本次执行绑定的来源/许可记录 |
| rho100 cell | 阻塞 | 计划要求使用经 C 核验的 `rho_ref(1250 K)`；D 当前可见登记表仍注明精确表/方程未独立取得，设计文件要求 C 的来源、单位、不确定度和 SHA 同步后才能确定 cell。故无法核实数值密度或 cell 参数 |
| 结构初始化 | 阻塞 | 方案只提出随机顺序放置；`r_min`、坐标生成器版本、坐标单位和 RNG/master-seed 派生尚未冻结 |
| MD 科学参数 | 阻塞 | 1600 K 预熔、降至 1250 K、NVT 及 0.5 fs 均在设计中标为提案或条件计划；没有可核验的正式 sampler/config 对象固定 timestep、温控器、参数、seed 与运行路径 |
| runtime | 部分可核 | `/workspace/onboarding/runtime.sh` SHA-256 `d8b6d0c05020eb9baf270ea17ef67c0f8c043d5b6d5f5eac692eba007969e926`。脚本配置 GPAW/MACE 环境，不含本次 Ag sampler 的启动配置；未据此运行计算 |
| 可见资源 | 只读已检查 | Linux x86_64；CPU 配额 4 个 CPU（cgroup `400000/100000`），可见 5 个逻辑 CPU；内存上限 32 GiB，检查时已用约 7.31 GiB；`/workspace` 可用 28,786,786,304 字节。无本任务运行中的 MD 进程 |
| 预检门禁 | 未通过 | 密度、生成器和运行参数均未冻结，sampler/许可证据不完整；按指令停止 |

## 试运行结果

- `ACTUAL_STEPS=0`；`SIMULATED_TIME_PS=0`。
- 初始坐标、位置/速度、RNG 状态、温度/能量记录、trajectory、checkpoint、运行 manifest 均未生成。
- 温度均值/波动、能量稳定性、RDF、Q6、配位数、MSD、异常短距离与原子丢失：未评估，因为没有轨迹或坐标。
- 20 ps 成本估算：无法从本次实测获得。积分步长/正式采样器仍未冻结，且没有首条轨迹的速度基准；不使用未验证假设外推。
- 本报告不把任何无序结构称为平衡液体。Ag64 step8000 的 `ARCHIVED_HOLD` 状态未改；其他五条轨迹未启动；没有 DFT、训练或角色分配。

## 解除阻塞所需证据

1. 可读的 C 密度核验记录及其来源、方程/表号、单位、不确定度、SHA，并明确 rho100 cell。
2. 冻结的生成器、算法版本、`r_min`、seed/RNG 派生和坐标约定。
3. 可执行且已审定的 sampler 及版本、势文件来源/许可、固定系综/温度调度/时间步长/温控参数、完整配置 SHA。
4. 与上述配置相符的实时 preflight 全项 PASS。

本报告只记录本次 preflight 阻塞和已完成的设计归档，不构成对另外五条轨迹、完整采样、DFT 或训练的授权。
