# 第 323 发失败现场

- 从第 322 发完成态启动：约 75.2 ps，4,627 个原子。
- 失败输入把新 Ag 放在 `(x,y,z)=(13.5159343, 24.0979196, 60.0) Å`。
- 与现有 Ag，ID 4598，坐标 `(13.5983764, 24.1212388, 59.5187412) Å` 的初始距离为 `0.489 Å`。
- LAMMPS 报告 EAM per-atom density 超表外推，之后在第 1,400 步（0.07 ps）报 `Lost atoms: original 4628 current 4627`。该尝试被终止。
- 失败的 input、console、log、dump、原始位置表以及发射前 data/restart 都保存在本目录。
- 失败前 `structures/state_current.data` 的 SHA-256 为 `4ccf887f5651370da9a0f319701d0b6d1313f5cceff2d1e0386a3f5c0a54ab83`；restart 的 SHA-256 为 `09bdefe422230863d2af0bcde2bc8e0a9793c399484e23800f05ae0edcf42d89`。

修正后从同一第 322 发状态续算：拒绝小于 2.5 Å 的入射初始近接，改为跨批连续的确定性随机流，并把非周期 z 边界设为 `m` 以保留溅射原子。参数修正不得被描述为原始 fixed-z 结果的无修改延续。
