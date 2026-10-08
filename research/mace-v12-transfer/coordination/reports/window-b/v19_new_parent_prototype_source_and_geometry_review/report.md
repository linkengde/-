# V19 新薄层原型的晶体来源与几何审查

## 结论

输入脚本、manifest 和三份输入 SHA256 均通过。三份几何均为52原子、PBC全开、无能量/受力标签；MAX 相基底上方有单层 Ag，基底上下终止与总结构最高原子层已分别辨明。现有公开训练36帧中找不到相同化学式结构；这支持几何/组成新颖性，但不证明真实表面稳定性或独立的物理源家族。

**尚不建议批准 DFT。** 任务卡要求确认空间群、晶格参数及原子位点的权威晶体学来源。仓库没有引用/来源记录；本次 Crossref 查询被网络代理 403 阻断，因此 a=3.07 Å、c=17.68 Å、Ti(4f) z=0.135、C(4f) z=0.567 仍是未核实假设。

## 分界面复算

|候选|基底下表面|基底上表面|Ag 最外层|Ag 到基底面间距|真空|最短 C–Ti|Ag 层面内最近距|
|---|---|---|---|---:|---:|---:|---:|
|AgC|Ti|C|Ag|2.40 Å|15.17 Å|2.132 Å|3.070 Å|
|AgSi|Ti|Si|Ag|2.60 Å|16.02 Å|2.132 Å|3.070 Å|
|AgTi|C|Ti|Ag|2.60 Å|15.18 Å|2.132 Å|3.070 Å|

所有 10 种元素对的周期最小距离都与 A 的 manifest 复算一致（最大差异见 audit.json）；含晶胞边界的距离已用 MIC 检查。Ag 层与所假定 a=3.07 Å 的基底恰好同格，几何上施加0%应变；这是构造设定，不是与实验 Ag 晶格的物理相干应变估计。

manifest 字段 `top_surface_species` 指 Ag 覆盖之前的 MAX 上表面终止元素；完整结构最高层实际都是 Ag。报告应使用“基底上表面终止”称呼 C、Si、Ti，避免把它误称整个界面的暴露终止面。单侧 Ag 覆盖和 Ti/C 不同的底表面终止构成非对称薄片，周期像会产生表面偶极耦合。几何尚未弛豫，不能从接触距离推断稳定吸附或平衡界面。

## 晶体学与独立性

生成器给出 P6₃/mmc (#194) 的常规12原子 MAX 相晶胞：Ti 2a + Ti 4f、Si 2b、C 4f；2×2 面内复制得到 C16Si8Ti24，再加四 Ag 得52原子。计量和 Wyckoff 轨道数自洽；但轨道坐标和晶格常数是显式写入的假设，生成器的自我复现不是外部结构证据。source_search 记录了未成功的检索，不将任何猜测文献列为已核实引用。

与训练36的组成分别为 C13Ag17Si8Ti13、C13Ag2Si10Ti23、C11Ag28SiTi9、C8Ag10Si2Ti12；三新原型均为 C16Ag4Si8Ti24，故不存在同式 exact/near 几何对。它们与旧小簇数据有明显的周期薄层拓扑差别。三者仍共享同一假定 MAX 相来源，只改变基底上终止与 Ag 接触类型；不能将三种候选夸大为三个独立母相或三种平衡形貌。

## DFT 前需完成

1. 由 A 提供可核验的 Ti₃SiC₂ 晶体学原始文献/数据库记录，核对 P6₃/mmc、a/c、Wyckoff 坐标及占位；若参数不同，A 在自己的输入目录修正并重新哈希，B 不改输入。
2. 说明所选 Ti₃SiC₂ 表面终止与单侧 Ag 覆盖的目标；报告一层厚度的有限尺寸、非对称薄片、固定未弛豫原子及偶极风险。几何不是平衡表面结构。
3. 先冻结 DFT 数值方案：PBE/PW500 eV；Fermi 0.1 eV；native extrapolated energy 与 free energy 分开保存；固定离子。2×2 面内、约15–16 Å 真空，建议在标签前进行小型 k 网格收敛（至少比较 Gamma 与2×2×1；若四 MPI 资源允许增加面内网格）；并比较偶极修正/真空收敛。不要将旧 Gamma 训练标签直接当作足够收敛证据。
4. 对照同一构型不同 k 网格的能量、力、界面力投影；确认四进程 MPI/ScaLAPACK、SCF、原子映射及 archive hash。A 需审查几何与参数来源后再冻结开发角色并登记 owner。

审查范围没有读取任何封存结构或标签，也没有启动 DFT、推理、训练、MD 或 TTM。几何筛查 PASS 只说明文件结构和周期几何自洽；晶体参数和物理结构资格仍待 A 解决。

Source/search evidence: {"attempted": "Crossref REST query for Ti3SiC2 crystal structure / lattice constants", "endpoint": "https://api.crossref.org/works?query.title=Ti3SiC2%20crystal%20structure%20lattice%20parameters&rows=5", "result": "UNAVAILABLE: network tunnel returned HTTP proxy 403; no source document retrieved", "repository_search": "No committed citation/table for these values found in research/mace-v12-transfer", "consequence": "P63/mmc, a=3.07,c=17.68,z(Ti4f)=0.135,z(C4f)=0.567 remain unverified assumptions; DFT readiness cannot pass crystallographic-source gate"}