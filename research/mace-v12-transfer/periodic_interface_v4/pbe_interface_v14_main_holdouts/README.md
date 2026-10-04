# 窗口 A 的 v14 预留盲测

两个输入分别为 AgC 和 AgSi，已生成并记录哈希；均不进入 v14 训练。先完成并核验 `pbe_interface_v13_holdouts/` 的当前两个作业，再确认其进程退出，之后运行本目录 `run_parallel_labels.sh`。脚本核对 window-a owner，只写本目录与 A 任务文件并逐点回传。四核环境中不要与当前 DFT 或 MACE 训练同时运行。

软件、分工和同步说明见 `../../coordination/README.md` 与 `../../coordination/ITERATION_PROTOCOL.md`。新窗口 B 不执行本目录任务。
