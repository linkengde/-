# B 窗口本地数据补充归档（2026-10-15）

本次只检查、复制和核验现有文件，未启动 DFT、MD、推理或训练，未修改原始结果和既有归档。

- 75 份本地结果、summary/progress 与已跟踪文件 SHA256 完全相同，保留对应路径于 manifest.json，不重复复制。
- 补充 25 份运行结束后的原始 GPAW 日志，保存在 raw_logs/。原始日志可能包含既有紧凑归档之后才输出的计时信息；逐文件前缀检查结果见 manifest.json。补充归档不替换紧凑日志，不修改原有校验文件。
- 补充用户上传的原始 COD 9009647 CIF，保存在 source/。与仓库现有 CIF 的差异仅为注释中的 COD 内部来源 URL，原文件字节予以保留；不把它当成新的独立结构或验证集。
- 5 个 state.gpw（约 2.31–4.14 GB/个）按既有 .gitignore 和 KEY_DATA_ARCHIVE_V19.json 的紧凑归档规则保留本地，未上传 Git。它们是重启二进制，不是本次补充的紧凑能量/力结果。manifest.json 列出路径及大小，未声称这些文件已远端备份。

manifest.json 记录审计时的 origin/main 提交、源路径、字节数、SHA256、既有归档位置和新增归档位置。SHA256SUMS.txt 覆盖本补充目录中除其自身之外的所有文件。

核验命令（在本目录运行）：

```bash
sha256sum -c SHA256SUMS.txt
```

MACE 小界面对照的模型结果已存在于 coordination/reports/window-c/three_model_static_comparison_20261010/，不在此重复归档。该报告使用的六个界面点是开发诊断，不能更改其数据角色。
