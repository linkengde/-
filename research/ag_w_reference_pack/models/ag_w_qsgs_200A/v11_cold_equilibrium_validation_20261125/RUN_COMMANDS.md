# 复现 V11 冷态平衡诊断

以下命令从仓库根目录执行。输出全部写入 `/tmp/agw_v11_cold_reproduction`；不会覆盖仓库内已归档结果或 V11 生成数据。当前归档运行使用仓库根目录 HEAD `40d6b5676c447c204026346957fee05714561817`。

```bash
BASE=research/ag_w_reference_pack/models/ag_w_qsgs_200A
OUT="$BASE/v11_cold_equilibrium_validation_20261125"
RUN=/tmp/agw_v11_cold_reproduction
DATA="$BASE/v11_generated_data/agw_v11_gyroid.data.gz"
POTENTIAL=research/ag_w_reference_pack/potential_files/derived/Dai_2006_EFS_AgW.eam.fs
LMP=/workspace/.venvs/ag-qsgs/bin/lmp
MPIEXEC=/workspace/.venvs/ag-qsgs/bin/mpiexec.hydra
export PATH=/workspace/.venvs/ag-qsgs/bin:$PATH
export LD_LIBRARY_PATH=/workspace/.venvs/ag-qsgs/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
mkdir -p "$RUN"

# 1) 独立 FIRE 力弛豫，最多 2000 步。
"$MPIEXEC" -n 4 "$LMP" -log "$RUN/v11_static_force_relax.lammps.log" \
  -var data "$DATA" -var potential "$POTENTIAL" \
  -var relaxed "$RUN/v11_force_relaxed.data" \
  -var dump "$RUN/v11_force_relaxed.dump" \
  -in "$OUT/in.v11_static_force_relax"

# 2) 0.5 ps 的共同 300 K NVT，然后 0.05 ps NVE。
"$MPIEXEC" -n 4 "$LMP" -log "$RUN/v11_cold_staged.lammps.log" \
  -var data "$RUN/v11_force_relaxed.data" -var potential "$POTENTIAL" \
  -var metrics "$RUN/v11_cold_staged_timeseries.dat" \
  -var restart "$RUN/v11_cold_staged.restart" \
  -var dump "$RUN/v11_cold_staged.dump" \
  -in "$OUT/in.v11_cold_staged"

# 3) 两次 0.5 ps 共同 NVT + 0.05 ps NVE 续跑，始终使用相同固定 ID。
"$MPIEXEC" -n 4 "$LMP" -log "$RUN/v11_cold_continue.lammps.log" \
  -var restart_in "$RUN/v11_cold_staged.restart" \
  -var restart_out "$RUN/v11_cold_continue.restart" \
  -var fixed_ids "$OUT/v11_cold_fixed_ids.include" -var potential "$POTENTIAL" \
  -var metrics "$RUN/v11_cold_continue_timeseries.dat" \
  -var dump "$RUN/v11_cold_continue.dump" \
  -in "$OUT/in.v11_cold_continue"

"$MPIEXEC" -n 4 "$LMP" -log "$RUN/v11_cold_continue2.lammps.log" \
  -var restart_in "$RUN/v11_cold_continue.restart" \
  -var restart_out "$RUN/v11_cold_continue2.restart" \
  -var fixed_ids "$OUT/v11_cold_fixed_ids.include" -var potential "$POTENTIAL" \
  -var metrics "$RUN/v11_cold_continue2_timeseries.dat" \
  -var dump "$RUN/v11_cold_continue2.dump" \
  -in "$OUT/in.v11_cold_continue"

# 4) 从约 1.65 ps 状态做 0.05 ps 高精度 NVE 积分器检查。
"$MPIEXEC" -n 4 "$LMP" -log "$RUN/v11_nve_precision_check.lammps.log" \
  -var restart_in "$RUN/v11_cold_continue2.restart" \
  -var fixed_ids "$OUT/v11_cold_fixed_ids.include" -var potential "$POTENTIAL" \
  -var metrics "$RUN/v11_nve_precision_check.dat" \
  -var dump "$RUN/v11_nve_precision_check.dump" \
  -in "$OUT/in.v11_nve_precision_check"
```

### 续跑和解释注意事项

- 每个 restart 都从前一段 NVE 结束状态读取，不能从最初数据文件重启续跑。
- 固定 ID 文件由力弛豫候选的底部固定层生成；任何变化都应视为另一组计算。
- `fix ave/time` 温度和能量数据的分辨率有限；审阅温度趋势时同时查看 LAMMPS 原始日志。
- `analyze_v11_cold_validation.py` 总结本目录归档数据。新复现输出在临时路径，需先独立审阅，避免覆盖已归档 JSON 和轨迹。
- 本流程没有 TTM 或激光 source。冷态未通过前，不应把它接作熔化模拟初态。
