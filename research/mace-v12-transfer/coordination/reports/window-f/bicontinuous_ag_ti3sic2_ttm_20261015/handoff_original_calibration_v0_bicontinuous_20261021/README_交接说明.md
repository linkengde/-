# 37,711 原子双连续 Ag–Ti₃SiC₂ 模型交接包（原 calibration-v0 势）

此包供另一个窗口基于原 calibration-v0 势继续建立/运行 LAMMPS 模型。它不包含 MACE V17 的运行检查点；换势函数后不能把 MACE 轨迹的速度、随机数状态或 checkpoint 当作原势函数的连续续算状态。

## 模型

- `model/AG_TSC_BICONTINUOUS_37711.data`：LAMMPS data，37,711 原子，4种类型，X非周期、Y/Z周期，110 × 98.0851344959 × 100.9207486328 Å。
- `model/AG_TSC_BICONTINUOUS_37711.extxyz`：OVITO/ASE 可读结构。
- `model/PHASE_ATOM_GROUPS.csv`：相/元素原子分组。
- `model/GEOMETRY_AUDIT.json`：几何与网络连通审计。
- `model/BICONNECTIVITY_SLICES.png`：形貌切片预览。

V2 几何包含 11,524 Ag 与 26,187 Ti₃SiC₂骨架原子，Ag质量分数约59.08%。这是基于独立理想晶格和随机场重建的三维双连续几何，非历史熔池轨迹复制。几何审计中，Ag最大连通分量约占Ag的99.94%；Ti₃SiC₂最大分量约占骨架的99.86%（2.75 Å几何邻接阈值）。这些只是几何审计，不代表势函数下稳定。

## 势函数

`potentials/potential.interface_calibration_v0.inc` 原样保留归档映射：Ag–Ag使用EAM，Ti–Si–C使用TSCBOP Tersoff，Ag–Ti/Si/C使用三项Morse交叉相互作用。元素/类型映射为 1 Ti、2 Si、3 C、4 Ag。include文件按相对路径读取 `potentials/Ag_u3.eam` 和 `potentials/TSCBOP.tersoff`。

calibration-v0 是原有探索性势映射，尚非生产级验证势。历史39,486原子模型在此映射下的短程测试曾报告局部高力/短Ti–C接触；这不等同于本37,711原子几何的结果，但本包没有提供“已稳定”结论。开始热模拟前应在该模型副本上做简短数值安全检查。

## 继续方式

新窗口从 `model/AG_TSC_BICONTINUOUS_37711.data` 建立独立的 LAMMPS 输入，include `potentials/potential.interface_calibration_v0.inc`。不要读取或叠加MACE V17权重，不要把MACE V17 checkpoint用于本势函数。请在新分支/新目录写输出，保留本包输入不变。

SHA-256见 `SHA256SUMS.txt`。本次只打包数据和势文件，没有在这里启动新的LAMMPS计算。
