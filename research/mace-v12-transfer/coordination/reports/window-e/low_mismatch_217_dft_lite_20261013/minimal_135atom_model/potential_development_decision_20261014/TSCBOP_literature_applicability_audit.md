# TSCBOP 文献与 Ti₃SiC₂ 适用性审计

## 结论与证据边界

**当前文件的物理适用性结论：FAIL / 不应作为已验证的 Ti₃SiC₂ 生产骨架势。** 文件内标题明确声称面向 Ti₃SiC₂ MAX 相；但原始论文及其晶格、键长、表面和应变验证范围无法在当前仓库和网络环境中独立核验。不能把文件标题或 LAMMPS 成功解析等同于物理准确。

| 项目 | 状态 | 依据 |
|---|---|---|
| 参数文件声称的材料 | `FILE_HEADER_MATCH` | 首行：“Bond Order Potential for Ti3SiC2 MAX Phase”；并列作者 G. Plummer、G.J. Tucker |
| 原始论文完整引文/DOI | `UNVERIFIED` | 仓库没有可核验的论文副本或 DOI。Crossref、OpenAlex 与 LAMMPS 文档访问均被网络代理拒绝（HTTP 403）；未推测或补写引文 |
| LAMMPS 文件格式/调用 | `IMPLEMENTATION_PASS` | 30 行三体元素参数表；历史静态 run0 使用 `pair_style tersoff` 并按 Ti/Si/C/NULL 映射。此项仅确认调用，不证明参数物理精度 |
| 元素、单位与截断映射 | `IMPLEMENTATION_PASS; ORIGINAL-CONVENTION UNVERIFIED` | include 使用 `units metal`、类型 1=Ti/2=Si/3=C/4=Ag；最大外 cutoff R+D=3.8838 Å。文件自身没有独立单位说明；原论文所用单位约定未核实 |
| 晶格常数重现 | `UNVERIFIED` | 没有已归档的 TSCBOP 晶胞优化、平衡体积曲线或与实验/PBE 晶格常数的直接对比。固定体相基准中明显的力残差提示警惕，但不能替代晶格常数测试 |
| Ti–C / Ti–Si 局部结构 | `FORCE_COMPATIBILITY_FAIL` | 已有同一 48 原子体相结构的力比较：Ti 与 C RMSE 约 2.1–2.2 eV/Å。未找到逐键长度、键角或声子谱验证 |
| 自由表面 | `UNVERIFIED` | 有一份从 Ag 界面移除 Ag 的 114 原子冻结 slab DFT/TSC-only 对照；它不是独立构造且弛豫的干净表面能/重构验证。公开本地档案没有原论文的表面测试范围 |
| 不同应变 | `PARTIAL / INSUFFICIENT` | 已有同一 48 原子父结构的四个 ±0.02 Å Ti/C 位移点。体积 ±2% 与剪切 ±0.5% 仍只是输入提案，未在 DFT 中完成；不存在充分的应变势响应检验 |

## 参数文件与计算映射

所审文件 SHA-256：`729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`，位置为 `research/mace-v12-transfer/coordination/reports/window-b/current_heating_potential_handoff_20261009/potentials/TSCBOP.tersoff`。对应旧 heating candidate 的 handoff 明确说 Stage69 原始配置身份未确认；因此本报告只认定这个哈希对应的实际候选文件，不把它追溯称为某个已认证生产版本。

135 原子 run0 使用的 `pair_coeff * * tersoff potentials/TSCBOP.tersoff Ti Si C NULL` 将 Ag 排除在骨架项外；`hybrid` 路由没有对同一对原子重复叠加。独立 run0 审计验证了原子 ID、坐标、PBC、力单位与文件 SHA；最大 R+D cutoff 3.8838 Å。以上支持“实现及映射正确”，并不说明参数可以重现 DFT。

## 实测材料表现

1. 同一势文件在 48 原子 Ti₃SiC₂ 体相固定结构上，与 GPAW/PBE/PAW 力比较的总向量 RMSE 为 **1.9965 / 1.9872 eV/Å**（两种 k 网格）。Ti 为 **2.2214 / 2.2110**，C 为 **2.1346 / 2.1246**，Si 接近零。两个 k 网格重现相同量级，不能解释成单纯 k 点噪声。
2. 在 135 原子接口派生的 114 原子冻结基底上，TSCBOP-only 对 TSC-only DFT 力 RMSE 为 **2.5624 eV/Å**；Ti、Si、C 分别 **2.6242、2.6072、2.4433 eV/Å**。严格内部配位样核心也有 **1.6265 eV/Å** RMSE。
3. 当前档案没有 TSCBOP 晶胞优化值、Ti–C/Ti–Si 平衡键长曲线、完整弹性张量、声子、表面能/重构或覆盖多个应变的独立验证，因此这些性质仍为 `UNVERIFIED`。

## 原始论文核查阻塞

参数文件只留有题名式描述和作者名，没有年份、刊名、页码或 DOI。本轮尝试访问 Crossref、OpenAlex 和 LAMMPS 的 `pair_tersoff` 说明，但 Codex Cloud 网络隧道均返回 HTTP 403，未取得外部文献原文或网页缓存。故本报告不编造引用，也不将“面向 Ti₃SiC₂”扩写成“论文已验证晶格、键结构、表面或应变”。后续若 A 路线仍考虑重标定，应先补齐可审计的原始出版物、参数来源及原作者测试数据，再决定是修正参数化还是替换错误/不匹配的参数文件。
