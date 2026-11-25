# Ag–W V11 冷态平衡与界面耦合专项验证

本目录是基于 V11 的独立诊断记录，基线为提交 `40d6b5676c447c204026346957fee05714561817`。本次新增工作只写入此目录；没有修改 V11 几何生成数据，也没有改动 V5–V10。

## 结论

**冷态基线：FAIL（尚未平衡）。**

- Ag、W 初速度分布、原子类型/质量映射和温度统计方式均核对一致；全体移动原子采用同一个 300 K NVT，不存在分相恒温器。
- 将初始静态最小化改为限步的力收敛尝试后，FIRE 在 2000 步上限停止，移动原子 RMS 力降至 `0.009995 eV/Å`，但最大力仍为 `1.239 eV/Å`，因此最小化没有完全收敛。
- 从该候选继续进行了三段各 0.5 ps 的共同 300 K NVT，并在每段后做 0.05 ps NVE。末段 NVT 结束时总体/Ag/W 温度为 `284.71/278.10/295.78 K`；随后高精度 NVE 末端为 `285.70/277.97/298.71 K`，W−Ag 温差扩大到 `20.74 K`。
- NVE 总能量漂移仅 `−4.17×10⁻⁴ eV / 0.05 ps`，但势能下降 `83.84 eV`、动能增加 `83.84 eV`，末端移动原子最大力/RMS 力为 `5.00/0.825 eV/Å`。能量守恒良好不代表结构已经平衡。
- 全部 `462,180` 个原子均保留，Ag/W 质量比为 `0.999995`；没有原子越过材料顶面 `z=200 Å`。仍有局部短距离与力热点，不能将当前构型作为已建立的 300 K 基线。

本轮没有启动 TTM 或激光热源。已有 TTM 软件联调记录仍然有效，但不等于冷态结构通过，也不构成熔化/凝固结果。

## 目录内容

- `V11_COLD_EQUILIBRIUM_REPORT.md`：完整方法、结果和判定。
- `RUN_COMMANDS.md`：从固定基线复现实验的命令。
- `in.v11_static_force_relax`、`in.v11_cold_staged`、`in.v11_cold_continue`、`in.v11_nve_precision_check`：LAMMPS 最小化、NVT/NVE 与高精度 NVE 输入。
- `v11_cold_fixed_ids.include`：续跑时重复使用的同一固定原子 ID 集合。
- `v11_*_timeseries.dat`、`*.lammps.log`：温度、能量、力及完整 LAMMPS 日志。
- `*.dump.gz`、`*.data.gz`、`*.restart.gz`：压缩后的结构、轨迹和检查点。
- `analyze_v11_cold_validation.py`、`v11_cold_diagnostics_summary.json`：可重复的分析程序与机器可读摘要。
- `SHA256SUMS.txt`：本目录文件校验值。摘要 JSON 另记录基线几何与 EFS 势文件哈希。

详细复现流程见 `RUN_COMMANDS.md`。压缩 restart 需先解压到临时目录再传给 LAMMPS `read_restart`。
