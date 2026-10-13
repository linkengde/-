# GitHub upload parts: Stage67 step1650 input package

本地完整文件夹：
D:\DUsers\DUsersdocument\新建文件夹 (2)\stage67_step1650_liquid_ag_complete_input_20261013

GitHub 分支中 ZIP 被拆成 4 个连续分块。下载全部四个分块后，按编号顺序拼接，得到完整 ZIP。

用 Python 拼接：
python -c "from pathlib import Path; p=['stage67_step1650_liquid_ag_complete_input_20261013.zip.part01','stage67_step1650_liquid_ag_complete_input_20261013.zip.part02','stage67_step1650_liquid_ag_complete_input_20261013.zip.part03','stage67_step1650_liquid_ag_complete_input_20261013.zip.part04']; Path('stage67_step1650_liquid_ag_complete_input_20261013.zip').write_bytes(b''.join(Path(x).read_bytes() for x in p))"

ZIP 大小：12,754,674 bytes
ZIP SHA256：e120ec43e02a6c7a0e0f574ae59e332e2074cc9cca322b33537be632b04125b4

分块校验值：
- part01：4,194,304 bytes；SHA256 6c64f12cdcd1df9a784e7a0dbb8a6db7b25d4fda59a61e4b89a0798bde5c29a9
- part02：4,194,304 bytes；SHA256 0ac52db8c49dcb5c705959581b07dda0630d4b38b05c7e05f43a6a997f40a8cb
- part03：4,194,304 bytes；SHA256 fee708c059892004bdfd35e3211c3b793dff57c05dcf0373bacf205abe1ae370
- part04：171,762 bytes；SHA256 8447ed2cb0e7e0de3aa6cbf18edb8944cb7b9d8061dd0c0e935fda9e915c1d38

ZIP 内含 15 个文件。解压后可按 SHA256SUMS.txt 校验。README 已注明现有数据缺少逐原子 PTM/Q6 标签和全量 LiquidCandidate 表。