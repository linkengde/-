# 候选02原始混合势与DFT固定结构诊断

保持当前交接的EAM+Tersoff+Morse参数不变，LAMMPS22Jul2025 update4、metal单位、单线程run0，计算16个已核验参考结构。原子按生成的LAMMPS ID还原输入顺序，再对应原始ID；核对元素/坐标/晶胞/PBC，允许仅周期包装坐标，未最小化、加热或积分运动。源势SHA256记录在report.json。

六个局部正负扰动点全原子力向量RMSE为1.870–2.607eV/Å，标记分离力误差2.618–5.606eV/Å，均超过项目0.05/0.10门槛。Ag–Si对应局部能差：混合势0.321158eV，DFT原生0.0275865eV、自由能0.0276292eV。Ag–C与Ag–Ti同样存在响应差异，详见report.json。以上只证明这些结构上的旧混合势与DFT不一致，不证明实际加热过程中全部能量释放的来源或大小。

完整混合势误差包含相内EAM/Tersoff和交叉Morse，不能把所有误差归给Morse，也不能用它补偿相内势误差。下一步应进行同几何相内/交叉项分解、Ag及骨架稳定性控制，再决定联合拟合或统一多元素势路线。交接配置未确认属于Stage69。

16点中10点是同几何数值控制，6点复用训练来源，均不是独立通过证据。不同绝对能量零点未匹配，因此未打绝对能量门槛分数；报告只列原始能量及同来源正负扰动的能差。有限展宽DFT原生/自由能分列，力一致性问题保留。没有拟合、采用新参数、恢复TTM或开启MD。

复现：`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 LD_LIBRARY_PATH=/workspace/.venvs/lammps-cross/lib /workspace/.venvs/lammps-cross/bin/python analyze.py`。predictions.json保存完整原子力供复核，metrics.csv是逐结构误差。输出文件会被复现脚本重写，仅在独立副本复跑以保留发布快照。
