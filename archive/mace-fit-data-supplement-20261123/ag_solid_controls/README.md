# 纯 Ag 固相 DFT 训练候选

该文件由当前仓库中 `research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/` 的已完成结果整理而来。

`curate_controls.py` 对每个选中结果逐项验证其 SHA-256 清单、原始输入哈希、GPAW 收敛日志、SCF 摘要、有限力、验证报告及方法设置，再将 `PW_PBE_energy_eV` 和 `PW_PBE_forces` 映射到项目训练字段 `REF_energy`、`REF_forces`。未改动源结果文件。

整理了 7 帧：Ag32 基准、沿 x 的 ±0.01 Å 原子位移、±2% 体积变化及 ±0.5% 剪切。网格为基准 8×8×8、位移 6×6×6、体积/剪切 8×8×8。方法均为 PBE/PW500、0.10 eV 展宽、GPAW 26.7.0；源控制的验证报告均为 PASS。

这 7 帧属于同一个 CIF 衍生 Ag32 晶体结构族，只作为未分配的训练候选。不能把它们分散到训练/开发/测试，也不能据此宣称液态或凝固行为已经覆盖。原始纯相控制登记的目的仍是数值控制；项目级训练集合并和划分审核尚未执行。
