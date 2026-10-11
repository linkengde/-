# Ag–Si DEV v2 registry DFT screening 放行审查

审查对象为 E 提交 `10e7acc27d203d7377c13f5de367690237e500e1`。结论：`BLOCKED`；本轮没有启动 DFT。

结构身份审查通过。E 包中的 R1–R4 各为 20 原子，组成均为 Ag7C4Si3Ti6；晶胞完全相同，PBC 均为 `[true, true, false]`。四个 `structure.xyz` 与 C 来源提交 `9a5482d93a8756850ab457a605cf75dfa1f8a3d4` 中对应的 XYZ 字节完全一致。源文件 SHA、E 包 SHA、各 `input.py` 内固定的结构 SHA、`PARAMETER_LOCK.json` 中的 SHA 全部匹配。E 清单中列出的文件 SHA 检查均通过。

几何参数审查通过。四个任务都是固定几何单点；原子与晶胞弛豫均关闭，脚本在计算前校验结构身份并在计算后断言原子坐标、晶胞及 PBC 未改变。脚本要求保存总能、force-consistent free energy、逐原子力和完整 3×3 应力张量，并检查 SCF 收敛及输出有限性。应力按包含真空的 slab 超胞体积归一，当前明确禁止把它作为体相 MACE stress 标签。

放行仍为 `BLOCKED`，原因不在 R1–R4 几何：E 的参数锁定只准备了 6×6×1 和 8×8×1 两套网格，两者都标为尚未生产批准；执行脚本还需要一个明确的正 MPI rank 数。另一个数据治理阻塞是脚本结果 JSON 没有写入或校验 `REGISTRY_SCREENING` 角色。当前源家族来自 COD 9009647/Kisi 1998，与此前评分的 Ag–Si 几何属于同一来源家族，所以这些结果只能用于 registry screening，不能转成独立 DEV，也不能进入 TRAIN 或 TEST。

我在 ASE 3.29.0 环境逐一运行了四个 `input.py` 的默认非计算预检。均返回 `INPUT_VALIDATED_DFT_NOT_EXECUTED` 和 `DFT_NOT_RUN`。未传入 `--run`，未执行 GPAW。

机器可读审查决定见 [`AG_SI_REGISTRY_DFT_RELEASE.json`](AG_SI_REGISTRY_DFT_RELEASE.json)。需要先确定并批准 k 点网格及 MPI rank，并把 screening-only 角色写入 release 与结果清单，再重新审查放行。
