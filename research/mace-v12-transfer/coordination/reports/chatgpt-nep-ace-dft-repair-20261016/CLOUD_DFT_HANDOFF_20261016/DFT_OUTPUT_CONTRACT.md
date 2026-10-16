# 输出合同及禁止误报

- 所有输入：见INPUT_MANIFEST.json；全部无新的DFT能量/力。
- 每个完成的DFT单点应保存：INPUT_SHA256、PAW四元素文件SHA256、GPAW/ASE/mpi进程数、PBE/PW500/k点/Gamma-centering/Fermi0.10/自旋/PBC、总能eV、free energy eV、32×3 forces eV/Å、收敛迭代数、墙钟、峰值内存若可查、SCF日志。
- 原始数据由历史GitHub仓库提供；本次交接由ChatGPT聊天界面生成；云端跑出的新DFT结果要标注执行节点与独立提交，不得混淆。
- 已完成：原子建模10个、冻结32原子1个、几何回读校验、离线预检、提交。
- 未完成：新的SCF、k点收敛、Ag液体DFT、完整骨架高温DFT、NEP/ACE训练、正式MD。
- 初步执行限制：先允许一个32原子新\`1×4×2\`SCF；禁止P1/P2自发并发/扩展、禁止覆盖。
