# 运行环境、依赖与复现

## 现有单温短测入口

在仓库根目录执行：

```bash
cd models/ag_ti3sic2_qsgs_1to1_finer
python -m pip install -r requirements-ttm-reference.txt
python -m unittest -v test_thermophysical_model test_qsgs_enthalpy_solver test_qsgs_source_integration
mkdir -p runs/single_temperature_handoff
python phase_field_qsgs/qsgs_enthalpy_solver.py \
  --spacing 5A --boundary A --heat-ps 0.1 --cool-ps 0.1 \
  --sample-ps 0.01 \
  --output runs/single_temperature_handoff/caseA_5A_heat0p1_cool0p1.csv
```

把 `--spacing 5A` 改成 `10A` 可做快速粗网格筛查；`--boundary B` 是理想热接触边界对照。程序会写 CSV 和相邻 `.json` 汇总。该入口从 300 K 零焓独立起算，不读取也不转换任何双温 restart。

检查完整 CLI 选项：

```bash
python phase_field_qsgs/qsgs_enthalpy_solver.py --help
```

程序当前自动根据显式扩散稳定条件选步长；CLI 没有可指定 `dt`、采样场快照、横向 tile 数、底部热沉或边界热通量的选项。热源经投影网格中心取样并乘 Ag 体积分数。虽然 `qsgs_source_integration.py` 已实现原始二值 QSGS 几何的精确解析积分，且 `QSGSEnthalpyGrid` 有对应构造参数，但当前 `qsgs_enthalpy_solver.py` 命令行没有暴露该参数。开展正式空间对照前需要补上该入口并测试。

## 依赖与本次核验环境

单温核心依赖见 `requirements-ttm-reference.txt`：NumPy `>=1.24`、SciPy `>=1.12`。本次交接核验环境：Python 3.12.14、NumPy 2.3.5、SciPy 1.17.0。无需 LAMMPS、原子势函数或 PyAMG 即可运行单温原型。

本次容器可见 5 个 CPU、约 34.9 GB 系统内存。资源信息只描述本次核验容器；目标云端服务器的实际核数、可用内存和 BLAS 线程数需运行前重新检查。已归档的双温 3×3 原生图 50 fs 运行约需 1.19–1.84 GiB 峰值内存、单组约 1.3–2.1 小时（具体随图模型/预条件器而异）；这不是单温运行时间估计，也不能直接外推到 1 ps。

## 复现数值基准

```bash
cd models/ag_ti3sic2_qsgs_1to1_finer
python phase_field_qsgs/stefan_benchmark.py
python -m unittest discover -v
```

完整测试集耗时和输出可能随环境改变。此次交接只运行了三个相关测试模块的 14 项单元测试；没有执行新的长时熔化/凝固计算。

## 建议的独立单温计算顺序

1. 保存运行参数、环境版本和清洁的新输出目录。
2. 运行现有短时 Case A 与 Case B 原型，确认焓/能量账本；勿覆盖仓库中旧 CSV/JSON。
3. 接入原始几何精确源积分及真实 Ag 连通输运，重新通过源能量、均匀导热、Stefan、局部输运和空间/时间收敛检查。
4. 冻结热源和边界假设后再做单温加热；满足熔化后单独运行关热冷却。
5. 所有阶段独立记录源能、总焓变化、潜热份额、液相体积和残差。不能把双温运行的温度或能量当作单温结果。
