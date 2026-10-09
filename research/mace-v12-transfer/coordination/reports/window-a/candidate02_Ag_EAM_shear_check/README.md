# 候选02 Ag EAM 剪切固定结构检查

用原始Ag_u3.eam、LAMMPS run0、metal单位和相同Ag基准/正负剪切k8 DFT结构。位置不动，不作积分/最小化。原子按生成ID排序后比较。三个输入文件SHA和势SHA在report.json。

两剪切点EAM全原子力RMSE约2.63e-5eV/Å，最大误差4.28e-5eV/Å。相对基准的EAM剪切能增量0.07579meV/atom，DFT原生增量0.06798meV/atom，差0.00781meV/atom；正负对称。表明该势在这两个纯Ag微小剪切点局部响应接近DFT。它不检验Ag–基底交叉项、绝对能量零点、弹性常数或热稳定性，也不确认Stage69配置。未改参数。

独立LAMMPS环境复现：`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 LD_LIBRARY_PATH=/workspace/.venvs/lammps-cross/lib /workspace/.venvs/lammps-cross/bin/python check.py`。请在副本中运行以保留报告快照。
