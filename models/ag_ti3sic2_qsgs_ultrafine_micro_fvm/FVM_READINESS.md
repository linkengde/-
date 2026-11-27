# 单温焓法FVM接入检查

## 两种可用输入

### main分支投影NPZ

求解器实际读取并经本轮逐值校验的字段：

- `ag_fraction`、`skeleton_fraction`：形状150³的双精度体积分数场；本几何是二值场。
- `ag_mask`、`skeleton_mask`：二值掩膜，互为补集。
- `shape_xyz`、`size_A`、`spacing_A`。
- 附加元数据：`size_m`、`spacing_m`、`seed`、QSGS参数和原几何SHA-256。

每个单元Ag与骨架体积分数之和严格为1。投影输入为 `ultrafine_500nm/fvm_inputs/*_fvm.npz` 和 `micro_10um/fvm_inputs/*_fvm.npz`。

### study分支原生NPZ

study分支单温焓求解器接受 `QSGSEnthalpyGrid(None, boundary, native_geometry_npz=path)`，原生几何字段为` skeleton`和`size_A`，由求解器把Ag定义为骨架补集。该接口也完成了三种子、两尺度、A/B边界的零热源短步检查。当前工作区及检查的仓库分支中没有找到名为`prepare_native_enthalpy_input.py`的准备程序或同名预制输入文件；本项目直接生成求解器接口所需的NPZ。

## 单位和数组方向

求解器尺寸字段用Å，转换关系是`1 Å = 1e-10 m`。

| 模型 | `size_A` | `spacing_A` | `size_m` | `spacing_m` |
| --- | ---: | ---: | ---: | ---: |
| 超细RVE | 100,000 Å | 666.6667 Å | 10 μm | 66.6667 nm |
| 微米RVE | 2,000,000 Å | 13,333.3333 Å | 200 μm | 1.333333 μm |

VTI将150³个体素保存为VTK CellData：`material_id`（0=Ag，1=Ti₃SiC₂）、`ag_fraction`、`skeleton_fraction`；图像点维度151³，原点0，spacing用米。验证器用VTK重新读回数组并校验x/y/z存储次序和单元数，避免把体素尺寸写成点间距或把轴顺序弄反。

## 边界和物理限制

主分支Case A是绝热骨架极限（骨架不贡献热容和导热）；Case B为两相理想热接触和共享温度，界面热阻为零。本轮未设置有限界面热阻、底面理想热沉、Marangoni力、电弧压力或液态流动。

零热源检查：投影输入12/12通过，原生输入12/12通过，温度恒定300 K、液相率为0、能量残差为0。按用户后续指示，对seed 1105两种尺度都进行了10步、0.5 ps/步的局部热源测试；峰值温度分别为314.10 K和300.736 K，最大Ag液相率为0，能量残差小于`8.3e-25 J`。热源为合成数值冒烟输入，不是实测工况。完整热史、跨越相变潜热的验证、网格收敛和真实热边界尚未做。

**接入判断：输入格式和零源状态GO；短程局部加热路径已冒烟通过；完整加热/熔化/流动计算NO-GO，直到热边界、材料参数和收敛条件有可追溯依据。**
