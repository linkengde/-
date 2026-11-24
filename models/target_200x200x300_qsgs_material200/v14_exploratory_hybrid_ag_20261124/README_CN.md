# V14 探索性混合势归档

本目录保存 Ag/Ti₃SiC₂ 的简化混合势候选、校准数据、短时稳定性筛查和 1024 原子升温—冷却试算。总体判断和限制见 [V14_EXPLORATORY_HYBRID_CN.md](V14_EXPLORATORY_HYBRID_CN.md)。

## 归档中的模型文件

- `fit_v14_targeted_ag_force5/`：四元素 V14 POD 力权重敏感性拟合候选。小界面 NVE 失败，不建议拿来做 MD。
- `v14_stage2_p0_ag_dft_diagnostics/ideal_cross_morse_fit/`：以 5 个界面构型校准交叉 Morse 项、以 1 个独立构型留出的结果。
- 小模型探索混合势：仓库公共势文件 `models/ag_ti3sic2_qsgs_1to1_finer/potentials/Ag_u3.eam` 和 `TSCBOP.tersoff`，加上存于各 `in.*` 文件中的交叉 Morse/ZBL 参数。

## 环境与复现

所需软件：带 EAM、Tersoff、Morse、ZBL 及 hybrid/overlay 的 LAMMPS；Python ASE、NumPy、SciPy；OVITO Python 模块 3.17.2 用于 PTM 分析。该目录内的原始 DFT 标签和 fitted POD 文件已归档，不需要重新拟合即可复查已有结果。

默认路径与本云端一致。换环境时可设置：

- `LAMMPS_BIN`：LAMMPS 可执行程序路径。
- `V14_POT_DIR`：包含 `Ag_u3.eam`、`TSCBOP.tersoff` 的势文件目录。
- `V14_AG_TARGET_K`、`V14_AG_HEAT_RAMP_PS`、`V14_AG_HOLD_PS`、`V14_AG_COOL_PS`、`V14_AG_FINAL_HOLD_PS`：小模型 Ag 指定温度历程。
- `V14_RUN_TAG`：输出目录后缀；重跑时建议使用新后缀，保留归档轨迹。

例如，在仓库根目录下用新标签复跑长保温探索算例：

```bash
V14_POT_DIR=/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/potentials \
LAMMPS_BIN=/workspace/venvs/podfit/bin/lmp \
V14_AG_TARGET_K=2200 V14_AG_HEAT_RAMP_PS=0.7 V14_AG_HOLD_PS=5.0 \
V14_AG_COOL_PS=2.0 V14_AG_FINAL_HOLD_PS=1.0 V14_RUN_TAG=2200K_replay \
python models/target_200x200x300_qsgs_material200/v14_exploratory_hybrid_ag_20261124/run_hybrid_1024_ag_melt_cool.py
```

这会再次生成约 1024 原子的周期理想界面测试；Ag 使用直接指定 NVT 温度历程，骨架全程 300 K 恒温。此输入是势函数筛查，不是 TTM 激光模拟。

## 文件校验

执行 `sha256sum -c SHA256SUMS` 检查本目录文件完整性。哈希清单相对本目录生成。

## 458,248 原子目标模型冷态验证

目标完整模型的原始结构、QSGS 构建输入和复现文件现归档在 `model_source/`。`Ag_Ti3SiC2_bicontinuous_200x200x300.data.gz` 解压后可供 LAMMPS 输入使用。该结构盒高为 300 Å，但实际材料原子只到 z≈198.4 Å；200–300 Å 是观察真空区。

`hybrid_458k_cold_screen/` 保存了成熟 Ag EAM + Ti₃SiC₂ Tersoff + 理想界面 Morse/ZBL 交叉项在完整 458,248 原子模型上的静态、FIRE、NVE 和分层冷态试跑。结论见 [完整模型冷态报告](hybrid_458k_cold_screen/SUMMARY_CN.md)。

固定底部 0–20 Å、让 20–40 Å 恒温层及 40–200 Å 牛顿层的骨架原子继续振动，是保留骨架传热通道的合理边界处理；但该候选混合势尚未通过冷态筛查：0.1 ps 无热源测试出现 Ti–C 1.439 Å 短接触，因此不应由这次测试批准加热。最大力小于 500 eV/Å 不能单独当作稳定性证明。

完整模型复跑前，先解压模型和 FIRE 松弛结构：

```bash
gzip -k -d model_source/Ag_Ti3SiC2_bicontinuous_200x200x300.data.gz
gzip -k -d hybrid_458k_cold_screen/cases/fire_fixedbase_100/minimized.data.gz
```

输入路径按本归档目录安排为相对路径；运行 LAMMPS 时保持当前 case 子目录作为工作目录。冷态分层测试没有使用 TTM 热源；恒温层是模型中的 Ag/Ti/Si/C 原子，当前模型未包含显式 Al 原子热沉。
