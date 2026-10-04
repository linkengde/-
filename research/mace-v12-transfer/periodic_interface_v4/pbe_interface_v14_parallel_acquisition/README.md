# 窗口 B 的辅助 DFT 任务

先读 `../../coordination/START_WINDOW_B.md`，成功领取 window-b 后安装缺少的软件并运行 `run_parallel_labels.sh`。本目录四个输入已准备好，前三个是采集点，最后一个 AgTi 点是 v14 预留盲测。每个标签核验后自动回传。输入生成器同时负责准备 A/B 初始任务包，已有输入时不要重复运行它。

本目录独立于当前窗口正在计算的 v13 AgC/AgSi，两组任务不能互相替代或重复启动。
