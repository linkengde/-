# 本窗口本地历史数据补传（2026-10-15）

基线：`8167544b05b013ba5bd6a47ceff68da0ab37fc66`。这是本窗口源目录中的历史数据快照；仓库主线已有 V19 工作，本归档不替换当前数据集、模型或验收结论。

对本地科研文件与基线仓库已跟踪文件按 SHA256 比对。`file_inventory.csv` 记录每个源文件的原始绝对路径、字节数、SHA256、仓库已有相同内容的位置或 ZIP 中的位置。相同内容只在 ZIP 中保存一次。

- 唯一补传内容：1015 个文件，解压共 11850567 字节。
- 直接可读的 V9–V11 补传副本：81 个文件，见 `priority_v9_v10_v11/`；其内容同时保存在 ZIP 中并由清单关联。
- `compact_history.zip`：保留 `lammps_ttm/` 与 `user_uploads/` 来源层级，涵盖早期 TTM/几何审计、势函数排查、小界面 DFT、MACE 试验、参考能和原始上传的 Ce/DOS 表等历史紧凑文件。
- `excluded_large_binary_inventory.csv`：187 个大型或二进制文件清单。这些文件未逐一比较仓库内容，不得据此断言全部未上传。

历史材料包含失败、未收敛、未选中或后来被修正的试验记录。文件名中的 PASS、valid/test、版本名、calibration 等保留原文，本次没有重新认定其科学有效性。历史脚本中的绝对路径属于原机器，后续运行前应检查；归档脚本仅作证据，不会自动执行。

本次仅进行文件内容哈希比对、缺失文件复制和归档完整性检查；未重新计算 DFT、训练 MACE、运行 LAMMPS 或 MD，也没有将历史结构加入任何当前训练/验证/盲测划分。

最大单文件纳入限制为 10 MiB；大型 GPAW 波函数、二进制模型/优化器检查点、生产模型、重启文件和大轨迹列在排除清单中，原件仍留在源目录。软件安装环境、编译树和缓存不属于本归档范围。

可从仓库根目录解压：

```bash
python3 -m zipfile -e research/mace-v12-transfer/coordination/reports/window-b/recovered_local_data_20261015/compact_history.zip /tmp/recovered_local_data_20261015
```

在本目录执行 `sha256sum -c SHA256SUMS.txt` 可核对上传归档及全部直接可读副本。ZIP 每个成员的 SHA256 记录在 `file_inventory.csv`。
