# E 节点势函数与 Codex 预检回执

采集时间：2026-10-11T03:16:08.444942+00:00。本回执仅上传势文件身份、解析摘要、节点资源与非科学性远端往返证据；不含 `Ag.eam.alloy`、PLT 参数表、PDF 或论文全文。

## 势文件与 ASE 解析

- 已解压本轮已上传 ZIP 到本机临时目录；未下载、替换或重新转换势文件。
- 本机文件：`/workspace/.ag-evidence/e-package.QlnLhK/Ag.eam.alloy`；大小 **780452 bytes**；SHA-256 `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`。
- 使用 C 分支 `f174cf6adbe35101f993a613f725586e27de8bc1` 的 `inspect_eam.py`（C 包 `SHA256SUMS.txt` 共 19 项全部通过），通过 E 的 runtime 环境 Python **3.12.14** / ASE **3.29.0** 解析为 **alloy**，元素 `Ag`。解析器完成有限性、元素元数据、网格和 cutoff 检查；没有创建原子、计算能量/力或运行 MD。
- cutoff 为 **5.995011000293 Å**，解析结果 `radial_tail_status=C_TAIL_REVIEW_REQUIRED`。两个登记候选盒子的半边长均小于 cutoff，解析器要求周期多重镜像，并标记径向尾部需 C 策略复核。此为解析诊断，不是采样放行。

## Codex 节点资源预检

- C `node_preflight_codex.py` 原始回执：`PREFLIGHT_STATUS=HOLD`，退出码 **2**；原因：`C_CONFIG_OR_OWNER_RELEASE_HOLD, REMOTE_SYNTHETIC_ROUNDTRIP_NOT_VERIFIED`。
- 本次资源子检查：4.0 CPU quota、5 个可见 affinity CPU；cgroup memory limit 34359738368 B、current 3859132416 B、headroom 30500605952 B；`/workspace` 可用 28795498496 B、inode 可用 2039326。Python 3.12.14 / ASE 3.29.0 / GPAW 26.7.0 / OpenMPI 5.0.7；runtime.sh SHA-256 `e32d19eef429939bcb3d0aa2d365d063b6dca16b4b4052b02da226fd7247e43b`。fingerprint SHA-256 `f6cd91beba645077ce89c5908fd24116053fa9266b06fe321ab537e52695e39c` 与当前 C HOLD 配置记录相符。
- 配置绑定的 E 采样根 `/workspace/.work/e/sampling` 尚未创建；持久化仍为 false。资源、身份、路径、势 SHA 子检查均通过。预检整体为 HOLD，因为 C 的配置仍未签署，且冻结配置里的 E synthetic roundtrip 字段仍为 HOLD；本轮不会修改 C 配置。

## 独立 GitHub 合成往返

- E 分支 `cloud-e-b1-dft` 的 fixture 上传提交：`649b7f4ef673ac618cf962f804f8081ff993ffc5`。远端 head 与本地该提交相符。
- 20,480 字节确定性合成 fixture SHA-256：`c3c582d73df54e9e9118efec0b96ca1c8faebba7b363ce3ebb007796cb250cf6`。在新的临时 Git 仓库执行 fresh `fetch`，再从远端恢复 fixture 和 manifest；恢复文件 SHA 与本地原件一致，字节比较通过。
- 仅 synthetic fixture、文件身份/解析摘要、资源预检和恢复收据进入此分支。`Ag.eam.alloy` 与参数/文献资料未加入、未上传。

```text
LOCAL_POTENTIAL_PASS=YES
ASE_EAM_PARSE_PASS=YES
NODE_PREFLIGHT_PASS=NO (C HOLD configuration; resource subchecks PASS)
REMOTE_ARCHIVE_PASS=YES (synthetic fixture only)
DFT/MD/TRAINING=NOT STARTED
```

等待 C 发布包含 E 当前 fixture 收据、节点新鲜资源和精确本机势 SHA 的新签署配置。所有科学计算继续 HOLD。
