# D Ag64 step 1000→2000 离线续算工具候选

状态：`PROPOSED_UNSIGNED_NO_MD`。本目录是执行门禁、checkpoint/归档协议和热化 QA 的离线设计与合成测试；不含 Ag 势文件，不接入科学采样器，不授权继续运行。

## 文件

- `REPORT.md`：结论与当前 HOLD 状态。
- `segment_guard.py`：独立、未接入生产执行器的门禁/原子归档候选。
- `test_segment_guard.py`、`SYNTHETIC_TEST_REPORT.json`：仅使用合成 checkpoint、功能模拟器、临时本地 bare Git 的测试。
- `PROPOSED_RELEASE_1000_2000.json`、`RELEASE_LOCK_UNSIGNED.json`：拟议身份与参数 pin。二者明确未签署。
- `SOURCE_ARCHIVE_VERIFICATION.json`：实际 D step-1000 归档的只读字节与数组核验结果。
- `THERMALIZATION_CRITERIA.md`、`RESOURCE_BUDGET.md`、`CODE_DIFF.md`：QA、资源及实现差异。
- `SHA256SUMS.txt`：交付文件校验。

源归档仍位于 D 的提交 `145ee8050cf95a3a1988e868c6b935e10538f723`；C 原始分析基准为 `d17f079d59268cf11071c9cbc7bd0ca08a26b8b3`。本目录不复制/改写 D 的 tar、checkpoint 或任何 A/B 文件。

测试命令：

```bash
python -m unittest -v test_segment_guard.py
```

完整测试通过只表示边界、身份、哈希、原子本地归档与本地 Git 读回的合成控制通过；不代表生产采样器已集成、D 端资源已实时确认、GitHub 远端恢复已执行或科学计算获批。
