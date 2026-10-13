# GitHub upload parts: Stage67 step1650 input package

The package folder is available locally at:
D:\\DUsers\\DUsersdocument\\新建文件夹 (2)\\stage67_step1650_liquid_ag_complete_input_20261013

Because the GitHub connector timed out on a single 12.75 MB archive, this branch stores the ZIP as four consecutive binary parts. Download all four files from this folder, concatenate them in numeric order, and verify the resulting ZIP.

Reassembly with Python:
python -c "from pathlib import Path; p=[''stage67_step1650_liquid_ag_complete_input_20261013.zip.part01'',''stage67_step1650_liquid_ag_complete_input_20261013.zip.part02'',''stage67_step1650_liquid_ag_complete_input_20261013.zip.part03'',''stage67_step1650_liquid_ag_complete_input_20261013.zip.part04'']; Path(''stage67_step1650_liquid_ag_complete_input_20261013.zip'').write_bytes(b''''.join(Path(x).read_bytes() for x in p))"

Full ZIP size: 12,754,674 bytes
Full ZIP SHA256: e120ec43e02a6c7a0e0f574ae59e332e2074cc9cca322b33537be632b04125b4

Parts:
- part01 — 4,194,304 bytes — SHA256 6c64f12cdcd1df9a784e7a0dbb8a6db7b25d4fda59a61e4b89a0798bde5c29a9
- part02 — 4,194,304 bytes — SHA256 0ac52db8c49dcb5c705959581b07dda0630d4b38b05c7e05f43a6a997f40a8cb
- part03 — 4,194,304 bytes — SHA256 fee708c059892004bdfd35e3211c3b793dff57c05dcf0373bacf205abe1ae370
- part04 — 171,762 bytes — SHA256 8447ed2cb0e7e0de3aa6cbf18edb8944cb7b9d8061dd0c0e935fda9e915c1d38

The ZIP contains 15 files. After extraction, verify the package contents with SHA256SUMS.txt. The package README notes that per-atom PTM/Q6 labels and a full LiquidCandidate table were not present in the available source package.
