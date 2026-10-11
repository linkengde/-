# 测试与静态核验结果

日期：2026-10-11 UTC。全部测试使用本地合成数据；没有执行材料采样或科学计算。

| 检查 | 结果 | 范围 |
|---|---|---|
| `python -m unittest -v test_codex_release_gate.py` | PASS，17/17 | 本地 bare-Git synthetic archive fresh-fetch、SHA 和恢复；签署/配置/代码哈希拒绝；owner/task 范围；远端通道、E fixture、势文件 receipt、存储模式及旧持久模式拒绝条件 |
| `python -m py_compile run_sampling_codex.py codex_release_gate.py archive_checkpoint.py node_preflight_codex.py inspect_eam.py` | PASS | 语法编译 |
| `python run_sampling_codex.py --help` | PASS | CLI 参数检查；未传入配置并未执行采样 |
| 任务清单对比 | PASS | 与授权执行包 `2959fd3e2396a7319bd1f5c7f4fec27ce442abe0` 的 6 项矩阵字节一致，SHA-256 `3232b70dde2c2c721b32114a9c30df280b6d0b64fd7c49d1aaca087fd0accd0b` |
| 旧 runner 对比 | PASS | `run_sampling.py` 与批准源提交中的原始字节一致，SHA-256 `b78b734798c9c075826e6611e3d2f7f156f0b05ff7c87159812ac76ed0e918a8`；`run_sampling_persistent_reference.py` 也相同 |
| 新 runner SHA pin | PASS | 代码实际 SHA 与 HOLD 配置 `execution_code_sha256` 一致；配置 SHA-256 `207451d2d43279b9a1ac2d661bbca56b60641890781090e9ea14ecb55ebbd34a` |

测试不验证节点当前状态、科学势函数是否合适、液态质量或长期存储服务，也不替代各 owner 的实时预检。E 的 owner fixture gate 仍是 HOLD；D/E/F 的本机精确 `Ag.eam.alloy` 收据仍缺失。
