# Stage67 原始模型生成链搜索报告

日期：2026-10-13  
扫描目录：D:\DUsers\DUsersdocument\新建文件夹 (2)  
范围：只读搜索文件名、相关文本和小型数据文件头；未运行 LAMMPS/MD，未修改原始文件。

## 结论摘要

**Stage67 的模型副本、区域文件副本、父模型基准和部分生成脚本均已找到。**

关键辨别：

- Stage67 清单里的 stage63_equilibrium_100.data 文件头显示 **150.1×200×200 Å、303,372 atoms、timestep=100**。后缀 _100 表示步数，不是 100×100 Å 尺寸。
- 150.1×100×100 Å 父模型是 source_baseline_150x100x100_300K.data，文件头显示 **75,843 atoms、timestep=1000**。
- 现存父模型基准可以直接作为 150.1×100×100 Å 冷态几何使用，也有几何复制输入及其验证报告。因此若“重建父模型”指从现有基准恢复/复制该父几何，**可以**。
- 若要求从更早 Stage61A 原始坐标端到端重新生成该父模型，**当前目录内证据不完整**：Stage63 脚本引用的原始 Stage61A 文件未按其原名找到。
- Stage63 重建链生成的是 150.1×200×200 Å 模型；它的 Stage61A 原始输入和中间 Stage63 data 当前未找到，无法仅凭这里的文件端到端复现历史 200×200 平衡快照。
- 已有 150.1×300×300 Å 的几何构建输入以 100×100 父模型为源，执行 Y/Z 方向 1×3×3 整数复制；这只生成几何，不代表加热或熔池验证。

## 关键文件与 SHA256

| 用途 | 文件路径 | 大小（bytes） | SHA256 | 备注 |
|---|---|---:|---|---|
| Stage67 使用的 200×200 模型快照 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage67_pressure_1fs_test_20260925\PRESSURE_P033\inputs\stage63_equilibrium_100.data | 38,909,611 | e193a7cde41a06fa04fd32ea9a80c60c4d5089d2a036a0c7fe6957ca57b4813f | 文件头：150.1×200×200 Å，303,372 atoms，step 100 |
| CONTROL 模型副本 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage67_pressure_1fs_test_20260925\CONTROL_P0\inputs\stage63_equilibrium_100.data | 38,909,611 | e193a7cde41a06fa04fd32ea9a80c60c4d5089d2a036a0c7fe6957ca57b4813f | 与压力支副本哈希相同 |
| Stage67 使用的区域文件副本 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage67_pressure_1fs_test_20260925\PRESSURE_P033\inputs\in.regions_x_surface.inc | 1,098 | e38b965a64cf2013720b3401f030ddab6fcfe57ebb5210b6373c85128ac8fb8a | 与 Stage67 输入清单哈希相同，定义 x 向分层和动态组 |
| Stage67 输入清单 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage67_pressure_1fs_test_20260925\PRESSURE_P033\inputs_manifest.csv | 4,334 | 8272d73657da8780e6d9dfc33696e88226355fa32753779c3ccb78a9fd54da1a | 记录来源路径、复制路径、SHA256 |
| 根目录同名区域文件 | D:\DUsers\DUsersdocument\新建文件夹 (2)\in.regions_x_surface.inc | 889 | 2199ac1ed0f0549e109b7605ea5ecd1f294c4cb28e3618193590b9c82f1951b2 | 与 Stage67 文件哈希不同，不可按文件名替代 |
| 150.1×100×100 Å 父模型基准 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage67_300x300_cold_start_geometry_20261013\input\source_baseline_150x100x100_300K.data | 9,001,578 | 3ec9b392fbe15974dd7b358d9ae8f00783ffa9b352e514b32ab904c10701c131 | 文件头：75,843 atoms，150.1×100×100 Å，step 1000 |
| 300×300 几何生成输入 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage67_300x300_cold_start_geometry_20261013\build_300x300_geometry.in | 323 | 007eefeb392f71e85fdfb859112267b4acb3b2fb2e61dff54db5d42f19a7e6e2 | 从父模型沿 Y/Z 做 1×3×3 整数周期复制 |
| 300×300 几何验证报告 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage67_300x300_cold_start_geometry_20261013\reports\GEOMETRY_VALIDATION.md | 2,380 | 708396f52de183b7ded82956f13b09dd77c86d1b59c9f736f534c68e53e91ed7 | 记录结构计数、盒子和几何检查 |

## Stage63 / Stage67 生成链脚本

| 环节 | 文件路径 | 大小（bytes） | SHA256 | 搜索到的内容 |
|---|---|---:|---|---|
| Stage63 几何准备 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage63_prepare_model.py | 4,697 | 0a5655377d4ea328aedc9e27c0e1693aa5221c8a63ace2e62d28ba130ed53451 | 读取 Stage61A 的 rotated_surface_x50_stage19_300K_passed.data，Y/Z 各复制 2 次，写 stage63_model_150x200x200.data |
| Stage63 输入文件生成器 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage63_make_input.py | 1,812 | b454e4287ee65030a9b7d15009b807ea941eced5c46bec3fbe4aad3e8985a8ff | 为 200×200 Stage63 模型配置平衡输入 |
| Stage63 100 步平衡输入 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage63_eq_reproduced.in | 1,032 | 018bb3764c8fecf1ace246ab4defa25a00fce73f638c5e9df75f4653d36a3e21 | 读取 200×200 Stage63 模型，执行 run 0 与 100 步，然后写 stage63_equilibrium_100.data |
| Stage63 TTM 输入 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage63_ttm_reproduced.in | 1,343 | 46f490507e488d74e852300af49c1e2cd671d09fe3b3fed2037409aef8a6dcad | 从平衡 restart 继续 TTM 500 步；本次未执行 |
| Stage63 启动脚本 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage63_run_reproduced.ps1 | 1,436 | 33a27c096a032e61d2c5690df08070f9b4c55162ec7211616b707eab4f464454 | 会先运行平衡输入，再运行 TTM 输入；本次未执行 |
| Stage63 OVITO 后处理脚本 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage63_postprocess_ovito.py | 10,100 | cba9737d18c66d96b92beefc2feb22da8676c3ba76c8002204ae883c9f184843 | 后处理辅助脚本 |
| Stage67 200×200 几何物化脚本 | D:\DUsers\DUsersdocument\新建文件夹 (2)\stage67_materialize_150x200x200.py | 6,301 | 754d5305b1c8b6997b03fa442a3c7b9f57d8e4791a5e59c798dc754b8e13199c | 另一个 100×100 到 200×200 周期平铺实现；脚本路径依赖桌面目录 |

## 缺项与可复现性判断

1. Stage67 输入清单记录模型来自外部路径 D:\AgTiProjectLink\Stage63_model_size_test\ttm_run_reproduced_equilibrium\stage63_equilibrium_100.data，哈希为 e193...813f；当前工作区中的 PRESSURE 与 CONTROL 副本均通过同一 SHA256 核对。
2. Stage67 输入清单记录区域文件来自外部路径 D:\AgTiProjectLink\Stage67_paper_arc_local_melting_100A\in.regions_x_surface.inc，工作区副本哈希与清单一致，为 e38b...fb8a。
3. 本次指定目录中未找到 Stage61A 原始结构文件 rotated_surface_x50_stage19_300K_passed.data，也未找到 Stage63 中间模型 stage63_model_150x200x200.data。Stage63 Python 脚本引用 D:\AgTiProjectLink 下的原始依赖；因此不能声称历史 Stage63 200×200 equilibrium 链已在此目录内完整可复现。
4. 150.1×100×100 Å 父几何文件本身存在，且有 SHA256 清单和已验证的 300×300 几何生成输入。可从此现成父模型进行几何复制；若要求从 Stage61A 的更早原始数据重新生成该父模型，则需要找回脚本引用的原始 Stage61A data 和配套模型来源记录。
5. 文件名 stage63_equilibrium_100.data 的实际盒子为 150.1×200×200 Å，不能把它误认为 100×100 父模型。
6. 根目录的 in.regions_x_surface.inc 与 Stage67 manifest 中的文件虽同名但哈希不同；应使用 manifest 对应的 Stage67 版本。

## 本次操作边界

只读扫描和读取文本/数据文件头，计算了上表文件 SHA256。未读取大型轨迹或 restart 内容，未运行 LAMMPS、MD 或 run 0，未改动原始文件。

