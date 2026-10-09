# 候选02 Ag EAM局部纯相检查

五个同方法k6已核验DFT控制，实际Ag_u3.eam、单线程LAMMPS run0，保持位置/晶胞/PBC不变，按生成ID排序原子力。原始势未修改。正负0.01Å位移的全原子力RMSE约0.0008884eV/Å，局部能差与DFT接近；这只支持该小位移控制，不是全面稳定性或独立界面验证。

体积+2%相对基准：EAM能差+2.13045meV/atom，DFT原生-6.66768、自由能-6.73348，趋势相反；-2%分别+2.30443、+11.5820、+11.6391。相对能量零点由同方法基准消去。体积k6尚未单独网格比较，故须先做k8数值对照，不直接判定最终误差或改EAM。原生/自由能保留有限展宽区别。未最小化、加热、训练或采用新势。

使用独立LAMMPS环境（LD_LIBRARY_PATH=/workspace/.venvs/lammps-cross/lib，OMP_NUM_THREADS=1，OPENBLAS_NUM_THREADS=1）运行check.py；在副本复现以保留发布输出。report.json含完整原子力及源文件hash。
