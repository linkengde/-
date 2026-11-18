# Round7 v11 Si 热点 DFT 标签归档

分支：`codex/round7-v11-si-hotspot-labels-20261118`。本分支从已归档的 v10 分支建立，保留 v10 归档，并新增 4 个已完成、已审查的 v11 Si 热点标签。

## 数据划分

| 索引 | 中心原子 ID | 角色 | SCF 迭代 |
|---:|---:|---|---:|
| 00 | 70588 | 训练候选 | 44 |
| 01 | 104359 | 开发验证留出 | 34 |
| 02 | 83874 | 训练候选 | 23 |
| 03 | 6223 | 训练候选 | 44 |

训练索引固定为 **00/02/03**；**01** 单独留作开发验证，不进入拟合。它不是盲测集。每个裁剪只用中心 Si 原子的 DFT 力作为目标；局部能量权重为 0，裁剪边缘力不作为 DFT 目标。

## DFT 与审查

方法为 GPAW PBE/PW350/Γ，Fermi–Dirac 展宽 0.10 eV；SCF 收敛条件为 density、eigenstates、energy 均为 `1e-5`。四个标签均达到三个收敛标记，坐标、原子顺序、中心 ID、有限力数组均通过审查。汇总见 `hotspot_label_audit.json`。

原标签程序在 DFT/extxyz 写出后，曾因将 NumPy 向量直接交给 JSON 编码器而没有生成部分摘要。归档的 JSON 摘要已从对应 extxyz 与 SCF 日志恢复；`elapsed_s` 留空，不伪造耗时。归档的 `label_hotspot.py` 已包含向量转普通列表的修正。

## 文件

`round7_pod_snap/models/v11_candidate/active_si_hotspots/` 包含四个标签的 extxyz、JSON 摘要、GPAW 日志、候选局部几何、提取清单和审查脚本。根目录的 `v12_label_split.json` 明确记录训练/验证划分与各文件 SHA-256。尚未审查通过的 v12 POD 拟合结果不在本归档中。

## 下载与校验

```bash
git clone --branch codex/round7-v11-si-hotspot-labels-20261118 https://github.com/linkengde/-.git
cd -/archive/round7-v11-si-hotspot-labels-20261118
sha256sum -c SHA256SUMS.txt
```

将其中的 `round7_pod_snap/models/v11_candidate/active_si_hotspots/` 复制到本机对应的 `round7_pod_snap/models/v11_candidate/` 下即可恢复标签包。
