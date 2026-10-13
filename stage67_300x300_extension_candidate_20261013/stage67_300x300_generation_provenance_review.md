# Stage67 300×300 Å 生成来源追溯与可行性审查

日期：2026-10-13  
范围：只追溯结构生成来源并判断同一流程能否支持 150.1×300×300 Å；未生成或修改模型文件，也未运行 LAMMPS、MPI、run 0、最小化或 MD。

## 结论

**能够确认 Stage67 step1650 的 200 Å 帧来自一个 Stage63 命名的 300 K data 输入；不能确认该输入如何生成，也不能确认 100 Å 到 200 Å 的具体结构变换。当前没有可复现的 Stage67 原始几何生成器、随机种子或变换脚本，因此不能按已验证的同一流程直接生成 300×300 Å 模型。**

本次没有生成几何候选。`stage63_equilibrium_100.data` 的文件名不能代替实际盒尺寸证据；该文件在当前工作区、已解包 Stage67 包及可见 Git 分支文件树中均未找到，无法用其清单 SHA 做实物校验。

## Stage67 来源链中实际找到的文件

审查基准为 handoff 分支 `cloud-d-stage67-300x300-extension-handoff-20261013`，commit `f4faca36e3987ddc577c49c6ce8c27a1731073da`；step1650 的已有审查分支起点为 `28b72301778b1f0cfc2bec11bb0148d69d743d2e`。

| 文件 / 记录 | 哈希或标识 | 能证明什么 / 不能证明什么 |
|---|---|---|
| `reports/Stage67_inputs_manifest.csv` | SHA256 `8272d73657da8780e6d9dfc33696e88226355fa32753779c3ccb78a9fd54da1a` | `model` 行记录源路径 `Stage63_model_size_test/ttm_run_reproduced_equilibrium/stage63_equilibrium_100.data`、期望 SHA256 `e193a7cde41a06fa04fd32ea9a80c60c4d5089d2a036a0c7fe6957ca57b4813f`，状态写作 `COPIED`。当前没有该 data 文件本体，故无法验证该 hash、元素数、盒尺寸或其生成过程。 |
| `input/implementation_reference/Stage67_200A_1fs_driver.py.reference.txt` | SHA256 `633a8a6a4d4927df8ee307f12febc4d8b8c18f6430b062d54fb43cd0507dbb2b` | 参考续跑 driver。硬编码 `BOX=(150.1,200,200)`、TTM=`15×20×20`、MD timestep=`0.1 fs`；在模型输入表中指向上述 Stage63 data，并在报告模板中称其为 200 Å盒的 300 K baseline。它记录的是既有 Stage67 连续运行/TTM 流程，不是原子结构生成器。 |
| `input/source_200A_step1650/high_frequency_step_1650.dump` | SHA256 `91d04a6fa3cb1b2fd9eb668b392cad04ecb593620723beac3e12ea4e7b8e85fc` | 已有演化帧：step 1650、0.165 ps、303,372 atoms、150.1×200×200 Å、`p p p`。它是运行输出，不包含初始生成脚本、随机种子或从前驱模型扩展的操作记录。 |
| `input/analysis/reference_manifest.json` | 包内已核验 | 将 step1650 与 Stage67 的 0.155–0.165 ps continuation checkpoint 关联；不提供初始结构生成方法。 |

在提供的 driver 中未找到 `replicate`、`create_atoms`、随机几何生成或结构相场初始化代码，也没有记录 Stage67 几何生成 seed。它使用绝对本机路径和既有 checkpoint，不能证明实际 100→200 Å 建模步骤。manifest 文件名含 `100`，但 driver 对加载后的 Stage67 data 声明盒长为 200 Å；由于源 data 缺失，无法判断 `100` 表示横向边长、阶段编号还是其他命名含义，也无法判断 200 Å 是源文件原尺寸还是由先前步骤扩展得到。

交接包的 README 与 `reports/missing_stage67_inputs.md` 同样将原始建模脚本、300 Å 结构规则列为缺项。包内可用的分析脚本用于 Q6/熔池判定，不生成 Stage67 原始结构。

## 相场 / 旋节分解代码：确实存在，但不是 Stage67 生成流程

在本地 Git 对象历史中找到独立的相场双连续形貌候选，提交 `17c038044a0b9c45da2c8f516c8d05a328527ad7`（2026-10-12，`Add phase-colored Ag-Ti3SiC2 deformation figure`）。该提交不属于当前可见的 `origin/main`、`origin/cloud-d-b1-dft` 或 Stage67 handoff/review 分支祖先；路径记录用于区分来源，不能当作 Stage67 文件。

| 文件 | SHA256 | 内容 |
|---|---|---|
| `outputs/AgTi3SiC2_phasefield_candidate_20261007/build_phasefield_candidate.py` | `dc59c2b62956e8fa0b90132f0b15efda670594a0f3ba5f87d75aec8e7b352971` | 守恒型 Cahn–Hilliard 相场和晶格位点掩膜映射前处理。 |
| 同目录 `audit.json` | `bfeee2d7af8952e634afd386d9ad0ad7861f655e821e2b5c08bb72dd6126dfb3` | 记录选定 seed/step、相区比例、形貌连通性及参数；明确声明是几何相场候选、未做能量/力计算或动力学。 |
| 同目录 `map_atomic_candidate.py` | `3dd64034dfde5d27c8ada711dcf9795c47b4534375d6ffb585842aa2a9a6b23f` | 按相场掩膜映射 Ti₃SiC₂/Ag 晶格位点并作局部硬距离投影；不是 Stage67 atom ID/坐标扩展脚本。 |
| 同目录 `inputs/planar_reference.data` | `80b53f26f0665607753ab7a7677d1d504e99b0d75ef78ec70f6d6c82cca7615a` | 上述相场候选绑定的 4,783 原子平面参考结构，不是 Stage67 step1650 或其 200 Å 前驱。 |

代码内可复现参数包括：网格 `36×36×48`；试验 seed `440101–440104`，审计记录选择 `440101`、step 500；任意形貌演化参数 `DT=0.002`、`MOBILITY=10.0`、`KAPPA_A2=0.45`；初始序参量均值 `0.20`、随机 Fourier 旋节不稳定波数带扰动标准差 `0.25`；TSC 目标掩膜体积分数 `0.60`；固体高度 50 Å、上方真空 20 Å。审计记录面内晶格矢量长度约 37.286 Å、夹角 120°，边界为面内周期、z 开放。相场自由能/迁移率/梯度系数没有针对 Ag–Ti₃SiC₂ 热力学标定。

因此它是**另建统计形貌候选**的代码依据，不是 Stage67 200 Å 模型的生成历史。当前参数和输入不支持直接设为 `X≈150.1 Å, Y=Z=300 Å, p p p`：面内尺寸由另一份 planar reference 推导，固体厚度与真空固定，z 为开放边界；其相区比例、晶格映射和坐标也不继承 Stage67 枝晶/骨架。将其修改成目标盒并不能称为原 Stage67 模型的精确扩展。

在 Stage67 包、当前工作区和可见分支文件树中没有找到“随机四参数法”的实际生成脚本、参数文件或 seed。该方法的参数目前均为 **UNVERIFIED / NOT FOUND**。上述相场代码内含随机 Fourier 初始扰动和四个 seed，但这不等于已找到名为“随机四参数法”的 Stage67 原方法。

## 两种建模目标必须区分

| 目标 | 当前是否可行 | 说明 |
|---|---|---|
| **原 Stage67 模型的精确扩展** | **否，证据不足** | 必须追溯实际前驱 data、它的盒尺寸和类型数、Stage63 生成脚本/参数/种子，以及从前驱到 200 Å 的确定操作。现有 200 Å step1650 帧是动态状态；即使找回一个 100 Å 初态并按整数倍复制得到 300 Å 初态，也不会因此得到与当前 0.165 ps 熔池帧同一条轨迹或同一形貌。 |
| **使用相同统计特征另建的新模型** | **方法上可另行研究；本次不能生成** | 可把独立相场候选作为研究起点，但须另行定义 300 Å 网格/边界、相区比例、晶格映射、密度和形貌指标，并校准统计参数。它不保留 Stage67 原区域坐标/原子 ID，也不能标为 Stage67 精确扩展。现有相场代码未达到该目标。 |

仅在拿到实际 `stage63_equilibrium_100.data`（SHA 与 manifest 一致）并证明原流程如何形成 Stage67 的 200 Å 初态后，才可以判断该前驱是否支持同一规则生成 300 Å **初始几何**。例如，若文件实测为 100×100 Å 且历史命令证实 Stage67 200 Å 正是从它作 `1×2×2` 周期复制，那么 `1×3×3` 可作为待审查的 300 Å 初始构型路线；目前这只是条件假设，不能执行或当作已验证结论。即使该路线成立，得到的也不是 step1650 熔池状态；生成对应 300 Å 熔池还需新的物理运行，超出本阶段授权。

## 若以后获批，最小复现验收项

1. 输入 data 与生成器、参数、seed 的 SHA 全部固定；证明原 200 Å 构型可以按记录流程复现。
2. 目标盒严格为约 `150.1×300×300 Å`、`p p p`；确认复制/边界操作与原流程一致，且不拉伸坐标。
3. 对精确扩展路线逐原子比较保留区 type、坐标和原 ID；记录新原子来源及新 ID 映射。
4. 从来源实测 Ti/Si/C/Ag 数量、质量/数密度和组成比例；与原模型按明确容差比较，不硬编码数量。
5. 检查 Y/Z 周期缝附近距离分布、跨缝邻接、空洞/重叠；与源模型内部邻距分布对照。
6. 对形貌报告枝晶/骨架的连通性、方向性、相关长度、截面积或体积分数。新统计模型只能报告与 Stage67 的统计相似度，不能声称逐原子同源。

本次只完成来源追溯。没有生成 `.data`、restart、LAMMPS input、TTM 映射或结构候选；没有修改 Stage67 文件、势函数、TTM、热源、边界或时间步。
