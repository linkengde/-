# Ti₃SiC₂ 四点骨架响应最终冻结模型审计

日期：2026-10-11（Asia/Shanghai）。B 数据基准：`34b8be988015af96d3fd8791fc041f6cc811e2b7`。

新增的 C 4f +0.02 Å 点补齐了碳原子的对称响应。原始 MACE-MP-0b3-medium 对 Ti/C 的局部 z 恢复响应分别低估约 20.6%/18.0%，V18 分别低估约 64.0%/72.0%。四点全原子受力误差仍显示 V18 的明显退化；这一证据继续支持 V19 从原始基础模型重新初始化，并对骨架能力设置回归门槛。此轮没有训练或修复模型。

## 四点共同基准与身份核验

读取了 B 的 `FINAL_RESULTS.md`、`final_results.json` 以及其中明确公开的四份归档。逐一重新运行现有 `verify_result.py` 的只读校验，并核对归档 `SHA256SUMS.txt`、B 最终 JSON 列出的各文件哈希及源输入哈希，全部 PASS。校验没有调用运行队列、发布器或 GPAW 计算器。

四点均有 48 个原子（24 Ti、8 Si、16 C），与父结构 `Ti3SiC2_baseline_k10x10x4_v19` 的元素顺序、显式 `local_new_id`、晶胞及 PBC TTT 完全相同；每点只移动指定原子的 z 坐标 ±0.02 Å。Ti 4f 的 ID 3 对应零基索引 2，C 4f 的 ID 9 对应零基索引 8。索引从归档 ID 数组和实际坐标变化共同解析，不将文件名中的 ID 直接当数组索引。

父结构源输入 SHA-256：`157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4`。父族沿用 `Ti3SiC2_bulk_parent_v19`，其文献结构来源为 COD 9009647；此前来源清单尚未确认远端 CIF 原始字节的独立下载哈希。本轮确认的是仓库既有父结构与四点归档的完整身份一致性，不把它提升为新的外部来源核验。

DFT 基准均为 GPAW 26.7.0、ASE 3.29.0、gpaw-data 1.2.1、PBE/PAW、PW 500 eV、电子展宽 0.10 eV、10×10×4 k 网格、SCF energy/density/eigenstates 1e-5、4 MPI ranks。方法 JSON 的规范化 SHA-256 均为 `702df71e586a61cd5785afc62653d7135a707032c91b72842b035b62198fd709`。原生 extrapolated `energy`（`force_consistent=False`）与 `free_energy` 分开保存并核对；本次响应比较直接使用同一归档的完整原子力。不同模型的绝对势能零点没有直接比较。

## 冻结模型与推理条件

两个序列化模型在读取前后均通过相同 SHA-256：

| 模型 | SHA-256 | 元素 / 截断 |
|---|---|---|
| MACE-MP-0b3-medium | `2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012` | 包含 Ti/Si/C/Ag；r_max = 6 Å |
| 冻结 V18 | `757586671174db85acbf7971ccfe92de19cc58378ee2431a797efa6d6163cb62` | 元素序列 C/Si/Ti/Ag；r_max = 6 Å |

Python 3.12.14、PyTorch 2.5.1+cpu、mace-torch 0.3.16、ASE 3.29.0、NumPy 1.26.4、SciPy 1.17.1，CPU float64，Torch/OMP/BLAS 2 线程、Torch interop 1 线程。显式使用 `eval()`；模型参数梯度关闭，力仍由坐标的自动微分求得。无优化器、权重更新或训练流程。

为统一四点条件，重新评估两模型×四结构共八次冻结推理。既有 Ti ± / C − 三点结果只用于交叉检查：DFT 力完全一致，模型力最大分量差小于 `4.8e-11 eV/Å`。新计算及归档只读校验合计约 7.65 秒；具体逐结构耗时记录于 CSV。没有重复任何 DFT。

## Ti/C 对称 z 恢复响应

采用 `Kz = −[Fz(+0.02 Å) − Fz(−0.02 Å)] / 0.04 Å`。相对误差为 `(K_model − K_DFT) / K_DFT`；负值表示偏软。

| 位点 | DFT Kz（eV/Å²） | 原始 MACE Kz | 相对误差 | V18 Kz | 相对误差 |
|---|---:|---:|---:|---:|---:|
| Ti 4f，ID 3 | 17.471248 | 13.868116 | −20.623% | 6.292036 | −63.986% |
| C 4f，ID 9 | 15.250488 | 12.511395 | −17.961% | 4.274081 | −71.974% |

| 位点 / 模型 | Fz(−0.02 Å)，eV/Å | Fz(+0.02 Å)，eV/Å |
|---|---:|---:|
| Ti / DFT | +0.308274 | −0.390575 |
| Ti / 原始 MACE | +0.180148 | −0.374577 |
| Ti / V18 | −4.166019 | −4.417700 |
| C / DFT | +0.393605 | −0.216415 |
| C / 原始 MACE | +0.373445 | −0.127011 |
| C / V18 | +2.792006 | +2.621043 |

V18 除中心差分响应偏软外，在两个 Ti 位移点均预测较大的负向力，在两个 C 位移点均预测较大的正向力，与参考力的符号变化不一致。仅报告斜率会掩盖这种静态力偏置，后续回归必须同时检查位点绝对力与对称响应。

## 全原子受力误差

力向量 RMSE 为 `sqrt(mean_atom(sum_xyz((F_model − F_DFT)^2)))`，最大误差为 `max_atom(norm_xyz(F_model − F_DFT))`，单位均为 eV/Å。该 RMSE 与把全部 3N 分量取平均的定义相差 √3，不能混用。

| 位移点 | 原始 MACE RMSE | 原始 MACE 最大误差 | V18 RMSE | V18 最大误差 |
|---|---:|---:|---:|---:|
| Ti −0.02 Å | 0.041095 | 0.128131 | 2.896862 | 4.474293 |
| Ti +0.02 Å | 0.039239 | 0.074294 | 2.880826 | 4.267399 |
| C −0.02 Å | 0.038849 | 0.060759 | 2.886734 | 4.279437 |
| C +0.02 Å | 0.040407 | 0.089406 | 2.891428 | 4.290535 |

四个等原子数构型的汇总向量 RMSE：原始 MACE 0.039908、V18 2.888969 eV/Å。逐 Ti/Si/C 元素及 x/y/z 方向误差见 `per_element.csv`，384 行逐原子结果见 `per_atom.csv`，8 行构型指标见 `per_configuration_metrics.csv`。这些结果属于同一父族的开发诊断，不能作为独立泛化分数。

## 对 V19 修复方案的新增证据

1. 原始基础模型仍为优先初始化候选；V18 仅保留为冻结退化对照。C 4f 新结果补强了此前 Ti 4f 的证据，但没有证明基础模型的骨架响应已达到最终精度要求。
2. 授权训练前，必须另行明确这些纯相数值控制是否、以及如何转为合法训练候选。当前角色严格保持 `numerical_pure_phase_control_not_training_or_independent_validation`；没有把四点复制进训练集或改变登记。
3. 受控修复的骨架回归应同时报告 Ti/C 的 Kz、位移两端绝对 Fz、全原子及分元素力误差，使用父族隔离的新开发结构决定训练停止条件。本次四点与基准共享父族且已经评分，后续不能重新标为独立测试。
4. 骨架回放、数据平衡与部分参数冻结仍是待授权对照策略。本次只增加诊断证据，未修改此前 V19 方案或报告，未宣布任何策略已经修复模型。
5. 此 Kz 是固定晶胞中特定位点的局部响应，不是宏观弹性模量，不代表完整 Ti₃SiC₂ 骨架已验证。亦不证明液态 Ag、熔化、凝固、界面高温响应或大规模 TTM 适用性。

## 交付与复现

- `response_comparison.csv`：Ti/C 两个完整响应对、模型力与相对误差。
- `per_configuration_metrics.csv`、`per_atom.csv`、`per_element.csv`：统一定义的受力误差。
- `source_identity.csv`、`source_verification.json`：逐输入/归档哈希、方法、父族、显式 ID、顺序、晶胞、PBC 及只读 verifier 结果。
- `inference_metadata.json`：精确软件版本、模型哈希/元素映射/截断、执行设置与耗时。
- `previous_result_crosscheck.json`：前三点与既有 C 输出的一致性检查。
- `reproduce_frozen_response.py`：仅访问明确公开四点的冻结推理脚本。指定新独立输出目录可复现，禁止指向 A/B 工作区。
- `validate_outputs.py`、`output_verification.json`：无需推理的输出身份、数量和指标核验。
- `SHA256SUMS.txt`：本目录其余交付文件的校验清单。

复现示例（仓库根目录；使用已经安装的 MACE 分析环境）：

```bash
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
PYTHONDONTWRITEBYTECODE=1 XDG_CACHE_HOME=/tmp/c-backbone-cache \
MPLCONFIGDIR=/tmp/c-backbone-mpl \
python -B research/mace-v12-transfer/coordination/reports/window-c/ti_c_backbone_response_final_audit_20261011/reproduce_frozen_response.py \
  --output /tmp/c-backbone-reproduction
```

本轮未读取封存结构/标签，未启动 DFT、MD、训练或 LAMMPS，未改动 A/B 工作区、任务、输入、checkpoint、进程或冻结模型。所有新增文件仅位于 C 独立报告目录；完成后等待用户授权。
