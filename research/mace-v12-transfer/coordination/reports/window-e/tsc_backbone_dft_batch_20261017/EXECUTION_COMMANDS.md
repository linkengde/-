# E：四个 Ti₃SiC₂ 表面固定几何 DFT 单点命令

执行锁 SHA-256：`678d1d287482e82039fc59e0c6961a92069a817f2c856c8745c82cf58bbd5863`

输入包/结构来自 E 提交 `cbc14c97fe5d962dc808a4cf48bca3e6b5588ea0`，C 审核来自 `669aafefd9efbda59b7ba088893e338d469d3f96`。每次只运行一个结构。执行器会重查资源、PAW/软件/结构 SHA、C 审核和分支状态；运行成功或失败后，只暂存该任务目录，推送 E 分支，再 fresh-fetch 并校验该目录的 SHA 清单。失败不会重试或自动改参数。

```bash
cd /workspace/-
git fetch --no-tags origin cloud-e-b1-dft
git switch cloud-e-b1-dft 2>/dev/null || git switch --track -c cloud-e-b1-dft origin/cloud-e-b1-dft
git pull --ff-only origin cloud-e-b1-dft

# Use the existing GPAW environment; do not install or upgrade packages here.
source /workspace/.venvs/gpaw-mpi/bin/activate
python -c "import gpaw, ase; print('GPAW', gpaw.__version__, 'ASE', ase.__version__)"

TASK_PACKAGE=research/mace-v12-transfer/coordination/reports/window-e/tsc_backbone_dft_batch_20261017
cd "$TASK_PACKAGE"
sha256sum -c SHA256SUMS.txt

# 先运行每项的只读预检。只有输出 PREFLIGHT_PASS 才运行下面对应的计算命令。
python execute_one_surface_task.py --task-id TSC100_TI0001_COD9009647_2x2x2_REF --lock-sha256 678d1d287482e82039fc59e0c6961a92069a817f2c856c8745c82cf58bbd5863 --confirm-5h-window --preflight-only
python execute_one_surface_task.py --task-id TSC100_TI0001_COD9009647_2x2x2_REF --lock-sha256 678d1d287482e82039fc59e0c6961a92069a817f2c856c8745c82cf58bbd5863 --confirm-5h-window

# 每项归档并 fresh-fetch 校验完成后，再对下一项重复这两条命令：
# TSC100_TI0001_COD9009647_2x2x2_RESP_Z020
# TSC100_SI0001_COD9009647_2x2x2_REF
# TSC100_SI0001_COD9009647_2x2x2_RESP_Z020
```

不要同时启动多个任务。每项固定使用 GPAW 26.7.0 / ASE 3.29.0、PBE/PAW、PW500 eV、FD 0.10 eV、标准 MP `(6,6,1)`、MSR1、默认 Davidson niter=2、SCF maxiter=160、MPI=4、OMP=1。固定坐标和晶胞；不做优化。输出的能量与应力只作诊断，获批用途是力标签候选。执行器没有 checkpoint；时间、内存或磁盘停止时保留日志并归档失败状态，不重启、不重试。
