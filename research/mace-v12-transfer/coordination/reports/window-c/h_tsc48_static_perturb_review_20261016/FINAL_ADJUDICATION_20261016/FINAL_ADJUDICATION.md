# H TSC48 单点候选：最终审查裁定

日期：2026-10-16（Asia/Shanghai）
被审查提交：`e923008453d71077421b651cf19c6705420d2356`
原始 C 审查记录：`../C_H_TSC48_SCIENTIFIC_REVIEW.md`

## 最终状态

```text
C_SCIENTIFIC_REVIEW = APPROVED
KGRID_DEFINITION = PASS
KGRID_CANDIDATE = 10x10x4
DATA_ROLE = HOLD
RESOURCE_CHECK = DEFER_TO_EXECUTION_WINDOW
RESOURCE_READY = NO
DFT_EXECUTION_READY = NO
execution_authorized = false
```

本文件记录用户给出的最终科学裁定，并补充原审查状态。原始审查报告及证据文件保持不变；本裁定不代表资源已通过，也不是 DFT 执行授权。

## 网格定义与适用范围

`KGRID_DEFINITION = PASS` 的定义为：候选的 `10×10×4` 使用与母体收敛对照相同的 GPAW 整数 Monkhorst–Pack 网格约定，不启用显式 Gamma override，也不增加 offset。它不是“显式 Γ 居中网格”的同义词。

母体的网格对照支持把 `10×10×4` 作为该候选单点的初始网格。扰动结构自身尚无 k 点收敛结果，因此本裁定不声称候选几何的网格收敛已被单独证明，也不扩展到应力收敛结论。

## 数据角色与资源边界

候选继续为 `DATA_ROLE = HOLD`。这次裁定不将它加入 TRAIN、DEV 或 TEST，也不改变同一 Ti₃SiC₂ 父族内既有数值控制数据的角色。若后续申请任何数据用途，仍需按父族整体审查。

资源审查推迟到执行窗口。执行前需核对当时的墙钟额度、内存与磁盘预算、checkpoint 写入和恢复、工作目录及归档通道。`RESOURCE_CHECK = DEFER_TO_EXECUTION_WINDOW` 不等于资源通过；当前 `RESOURCE_READY = NO`。

当前 `DFT_EXECUTION_READY = NO` 且 `execution_authorized = false`。本文件没有启动计算，也没有授予启动权限。唯一单点仍需用户另行明确授权，并在执行窗口通过资源检查后才能运行。

## 来源

- H preflight 报告、候选结构、母体来源、k 点比较和既有标签证据见原始审查目录的 `evidence/`。
- 用户最终裁定：2026-10-16 对话指令，内容如本文件“最终状态”所列。
- 未进行 DFT、MD、AIMD、TTM、松弛或训练。
