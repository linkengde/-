# D Ag64 暖机离线审计补充包

本目录是对 D 已完成 step-1000 暖机的**事后离线补充**。它不改变 D 的历史提交、原始 `warmup-*.tar.gz`、checkpoint、归档 manifest 或远端恢复回执；其中任何 README、图、诊断 CSV、分析代码和输入副本都不是原始 tar.gz 的成员。

## 权威 D 归档与文件对应

- D 分支：`cloud-d-b1-dft`
- D 的恢复校验提交：`145ee8050cf95a3a1988e868c6b935e10538f723`
- D 原归档提交：`11803909226538ab83847d1f06827a8292501139`
- D 原归档路径：`codex_warmup_archive/C-D-AGL-1250K-rho097-seed01/b62a818e58a2cf42/`
- 原 tar.gz：`warmup-b62a818e58a2cf42.tar.gz`，100700 bytes，SHA-256 `b62a818e58a2cf423e4cb1b9ed0468ae73ad578a75c4b2cb9fb3f28baf75e7b1`
- 原 `archive-manifest.json`：SHA-256 `ab9f0dffa8ef145eac7792dfe1e14e07c1354df1c2152dafec97bad1dc3a57b1`
- 原 `REMOTE_READBACK_RECEIPT.json`：SHA-256 `8bc2b17738f5e214098964c0635ef3505134726f8b4d8e7323cd4c5b03591bd4`

`archive-manifest.reference.json` 与 `remote-readback-receipt.reference.json` 是以上 D sidecar 的字节级副本，方便核对，不是 tar 内成员。`inputs/from_D_archive_members/` 是从 D 的 fresh-fetch 恢复目录逐字节复制的18个归档成员；`archive-manifest.reference.json` 验证这些成员的路径、文件长度和 SHA-256 全部通过。该副本没有包含 `Ag.eam.alloy` 或论文文件。原始归档仍以 D 分支的 tar.gz 为权威对象。

根目录的 `thermal.csv` 是原归档热日志的逐字节副本；其 SHA 与 manifest 中 `thermal.csv` 条目相同。新增热化统计、结构统计、图表和建议均由这些只读输入计算。参见 `REPORT.md` 的方法、定义及判定边界。

## 可复现检查

使用 Python 3.12、ASE 3.29、NumPy、SciPy 与 Matplotlib。分析程序不会创建 calculator，不运行动力学、DFT 或训练：

```bash
python analyze_archive.py \
  --archive-dir inputs/from_D_archive_members \
  --manifest archive-manifest.reference.json \
  --outdir .
```

分段门禁原型测试只用合成字典，不接触模拟器。它不是生产执行器，也不授予计算权限：

```bash
python -m unittest -v test_segment_gate_proposal.py
```

当前分段提案仍为 `PROPOSED_NOT_SIGNED`。即使门禁单元测试通过，也不能启动 MD；实际执行前必须由 C 对新执行器和新配置重新审查、分别固定 SHA-256 并正式签署。

## SHA-256 清单

`SHA256SUMS.txt` 覆盖本目录中除该清单本身以外的全部文件，路径相对于本目录。校验命令：

```bash
sha256sum -c SHA256SUMS.txt
```

清单自身由包含它的 Git commit 内容固定。D 原 tar 的 SHA-256 在上文和 `REPORT.md` 中单独记录，并由 D 的原 `archive-manifest.json` 与 fresh-fetch 回执支持。
