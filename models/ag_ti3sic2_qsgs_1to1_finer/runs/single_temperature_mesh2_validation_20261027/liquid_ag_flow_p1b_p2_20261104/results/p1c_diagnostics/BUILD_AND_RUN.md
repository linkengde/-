# P1c 复现说明

使用与 P1b 相同的 Basilisk 快照：comphy-lab/basilisk-C v2026-07-20，提交 9bb4f30aaf195feaa7a3a613468a1c888f5a4b9e。先设置 Basilisk 的 src 路径，例如：

    export BASILISK_SRC=/path/to/basilisk/src
    export BASILISK=$BASILISK_SRC
    export OPENGLIBS=-lfb_tiny

编译动态诊断和轴对称面力检查：

    $BASILISK_SRC/qcc -disable-dimensions -O2 -Wall -autolink p1c_vof_dynamic.c -o p1c_vof_dynamic -lm
    $BASILISK_SRC/qcc -disable-dimensions -O2 -Wall -autolink p1c_axisymmetric_force_balance.c -o p1c_axisymmetric_force_balance -lm

动态程序参数依次是：每液滴半径单元数、力项模式（0 联合、1 法向变量项、2 切向项）、AMR 开关（0/1）、毛细步长倍率、切向离散诊断系数。联合标准基线用系数 1.0；0.9792 只复现本报告中的静态力平衡候选修正，不是建议的生产参数。

    ./p1c_vof_dynamic 8 0 1 1 1.0
    ./p1c_vof_dynamic 16 0 1 1 1.0
    ./p1c_vof_dynamic 32 0 1 1 1.0
    ./p1c_vof_dynamic 16 0 1 0.5 1.0
    ./p1c_vof_dynamic 16 1 1 1 1.0
    ./p1c_vof_dynamic 16 2 1 1 1.0
    ./p1c_vof_dynamic 8 0 0 1 1.0
    ./p1c_vof_dynamic 8 0 1 1 0.9792
    ./p1c_vof_dynamic 16 0 1 1 0.9792
    ./p1c_vof_dynamic 32 0 1 1 0.9792

轴对称面力积分：

    ./p1c_axisymmetric_force_balance 8
    ./p1c_axisymmetric_force_balance 16
    ./p1c_axisymmetric_force_balance 32

原始逐步输出和 stderr 应保留；CSV 中的压力残差是 Basilisk mgp.resa 原始量，同时给出其相对本步压力 tolerance 的比例，避免跨不同时间步直接比较残差数值。
