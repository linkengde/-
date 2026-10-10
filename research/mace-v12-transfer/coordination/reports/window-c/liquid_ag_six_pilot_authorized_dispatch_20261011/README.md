# 六个液态 Ag pilot 执行包

C 分支：`cloud-c-liquid-ag-six-pilot-execution-20261011`。先验证 `SHA256SUMS.txt`，再读 `REPORT.md`、`C_RELEASE_GATES.md`，各机分别读 `START_D.md`、`START_E.md`、`START_F.md`。

用户已授权六条液态候选轨迹、必要的同几何 Gamma/2×2×2 数值对照及合格后的六个 GPAW 单点；D/E/F 各两项。G/H 待命，其余24项 HOLD。无需再次批准整轮，科学和资源门槛仍须通过。

密度论文与 WMH2006 论文／NIST 发布说明已核验；实际 `Ag.eam.alloy` 已收到，四参数文件SHA与下载记录一致，ASE3.29解析和近零尾段审核通过。`SOURCE_RELEASE.json` 已冻结来源；存储／暖机／DFT运行门仍关闭，不能直接启动采样或 DFT。

工具已备好：`preflight_owner.py` 只读采集本机证据；`inspect_eam.py` 解析实际势文件；`run_sampling.py` 分暖机／续算并保存重启状态；`qa_sampling.py` 和 `qa_independence.py` 审查真实轨迹；`run_dft_pilot.py` 在放行后启动本机四rank受监督单点；`compare_kpoints.py` 比较已完成的同几何结果。`C_RELEASE_CONFIG.json` 与 `DFT_RUN_TEMPLATE.json` 都是未放行模板，不是可运行科学输入。

所有代码验证仅限编译、帮助、纯函数数值检查及关闭门的拒绝执行检查；没有运行本六任务的 MD/DFT，没有候选坐标或 SCF 结果。真正运行须使用 C 随来源／资源审核发布的配置 SHA。参数表、完整轨迹和 checkpoint 存于各机独立持久根，不上传本仓库。
