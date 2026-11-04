# P1d CLSVOF + integral.h 低 Ca 基准复现

本目录使用和 P1b 相同的 Basilisk 镜像与依赖。源码版本、哈希和编译器信息见 p1d_source_manifest.json。设置 Basilisk src 路径后执行：

    export BASILISK_SRC=/path/to/basilisk/src
    export BASILISK=$BASILISK_SRC
    export OPENGLIBS=-lfb_tiny
    $BASILISK_SRC/qcc -O2 -Wall -autolink marangoni_Ca0066.c -o marangoni_Ca0066 -lm
    timeout 900 ./marangoni_Ca0066 > Ca0066_stdout.log 2> Ca0066_stderr.log

源码以归档的官方 src/test/marangoni.c 为基线，仅将 Ca=0.66 改为 Ca=0.066，并在 3*t0 测量后加入停止事件。官方测试内部运行 LEVEL=7、8、9，即 8、16、32 cells/radius。不要删除输出日志；stdout 中包含逐步速度、质量漂移和多重网格记录，stderr 输出每档最终速度对照。


运行从温度场读取 σ 的 R16 接口烟雾测试：

    $BASILISK_SRC/qcc -O2 -Wall -autolink results/marangoni_temperature_sigma_R16.c -o results/marangoni_temperature_sigma_R16 -lm
    timeout 180 results/marangoni_temperature_sigma_R16 > results/temperature_sigma_R16.stdout.log 2> results/temperature_sigma_R16.stderr.log

此测试中的 T 是解析指定的静态线性场，不是热传导求解结果。
