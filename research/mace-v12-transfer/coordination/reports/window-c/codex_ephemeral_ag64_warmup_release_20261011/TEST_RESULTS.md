# 测试与静态核验结果

日期：2026-10-11 UTC。所有本地集成测试使用合成数据；没有执行 Ag 原子采样或其他科学计算。

| 检查 | 结果 | 范围 |
|---|---|---|
| `python -m unittest -v test_codex_release_gate.py` | PASS，18/18 | 本地 bare-Git synthetic archive fresh-fetch、SHA/恢复；未签配置、错误 SHA、代码 pin、owner/task 错配、缺远端通道、E fixture 状态、势文件 receipt、存储模式及旧持久模式拒绝；另验证当前签署配置仅含 D 首项 1000 步、0.5 fs、E/F host release=false，且只对准确势 SHA 声明 C 已批准的窄尾部策略 |
| `python -m py_compile run_sampling_codex.py codex_release_gate.py archive_checkpoint.py node_preflight_codex.py inspect_eam.py` | PASS | 语法编译 |
| `python run_sampling_codex.py --help` | PASS | CLI 参数检查；未提供配置或启动采样 |
| 任务清单对比 | PASS | 与授权执行包 `2959fd3e2396a7319bd1f5c7f4fec27ce442abe0` 的六项矩阵字节一致，SHA-256 `3232b70dde2c2c721b32114a9c30df280b6d0b64fd7c49d1aaca087fd0accd0b` |
| 原 runner 对比 | PASS | `run_sampling.py` 及 `run_sampling_persistent_reference.py` 保留批准源提交字节，SHA-256 `b78b734798c9c075826e6611e3d2f7f156f0b05ff7c87159812ac76ed0e918a8` |
| 新执行代码 SHA pin | PASS | 五个执行文件实际 SHA 与当前配置 `execution_code_sha256` 相符：archive `dd4a99a7…eecfb2`、gate `82c3a088…b0c5e3`、inspect `d3a43758…526a77`、preflight `cbd10c7f…6c46b5`、runner `6b54d094…b943583` |
| 径向尾部审查 | PASS，参数元数据/样条检查 | 精确用户提供势 SHA `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`；ASE 3.29.0 的一个 `dr` 内插值最大绝对密度 `1.9791e-14`，最大绝对 pair phi `2.0657e-13 eV`。未计算原子能量/力或运行采样；不推广到其他文件/解析器 |
| 配置完整性 | PASS | 现行 C 补充签署配置 SHA-256 `2e122b258a48fe54bf603709b5505a301c1257f398c5625eecc2ca94998c0052`；前一签署配置 SHA-256 `4aee9d64859736f35459259625598f767d5f3e0d8b33bc71f0abb5baf31639e2`；最初 HOLD 配置 SHA-256 `207451d2d43279b9a1ac2d661bbca56b60641890781090e9ea14ecb55ebbd34a` |

节点证据的范围、owner 提交与逐文件 SHA 见 `OWNER_EVIDENCE_INDEX.csv`。测试不替代 D 当前会话用签署配置 SHA 的实时 node preflight；仅当其返回 PASS，D 才可运行被签署的唯一暖机 task。测试也不验证势函数对熔点、扩散或凝固物性的准确度。E/F 暂不放行。
