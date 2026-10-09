# 当前候选02加热实际势文件交接

供 window-a 的参考数据核对及后续旧势能量—力评估使用。已读取 origin/main 9d4de7c 的首轮报告。本包提供当前候选02加热实际使用的文件，**不能标称为已确认的 Stage69 势**。Stage69 几何报告提到 ZBL 就绪，但未附实际输入及势文件；本包无 ZBL。两者身份仍需历史原始文件证明。

## 可复现配置

工作目录为本目录；读取 potential.interface_calibration_v0.inc。units metal；atom_style atomic；类型 1=Ti、2=Si、3=C、4=Ag。实际 pair_style/pair_coeff 和三组 Morse 参数均完整保留在 include 文件，EAM/Tersoff 文件按原相对路径提供。Morse 截断 8 Å；include 中无 pair_modify shift yes。不要静默添加能量平移。源路径和文件 SHA256 见 manifest.json。

此配置为 calibration v0，未通过生产验证。当前热分支因交叉势能量释放异常暂停。本次没有拟合、没有换势、没有启动 MD，也未改变骨架、比例或 TTM。

## 历史间距数据核对

historical_classical_rigid_scan.csv 为上传包 Stage70 的经典混合势刚性扫描，18 组终止/配准、每组 27 点，共 486 点。字段明确是 Wad_rigid_model_J_m2，远距 d=12 Å 作分离参考。它不是 DFT 原始扫描，不能当独立真值拟合其自身。论文表格的 18 个 d0/Wad 锚点也不是逐距离能量/力数据。Stage70 重建几何未证明与论文原始 slab 相同，不能把旧误差表当严格同构型验证。

下一步可用本包计算当前候选02旧势对已核验16案例的能量/力响应，并明确自由能/外推能量口径。绝对黏附能仍需匹配的分离参考和界面面积定义。Stage69 专属误差计算仍待其真实文件，不用本包替代其身份。
