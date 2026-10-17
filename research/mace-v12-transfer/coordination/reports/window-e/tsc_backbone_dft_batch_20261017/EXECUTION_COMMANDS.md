# E：四个 Ti₃SiC₂ 表面固定几何 DFT 单点命令

执行锁 SHA-256：`579f5803cd6c3e5d522ecf396c4f25fb06ea9fe8939f9ef684e750702f9edb46`

输入包/结构来自 E 提交 `cbc14c97fe5d962dc808a4cf48bca3e6b5588ea0`，C 审核来自 `669aafefd9efbda59b7ba088893e338d469d3f96`。每次只运行一个结构。执行器会重查资源、PAW/软件/结构 SHA、C 审核和分支状态；运行成功或失败后，只暂存该任务目录，推送 E 分支，再 fresh-fetch 并校验该目录的 SHA 清单。失败不会重试或自动改参数。

```bash
cd /workspace/-
git fetch --no-tags origin +refs/heads/cloud-e-b1-dft:refs/remotes/origin/cloud-e-b1-dft
if git show-ref --verify --quiet refs/heads/cloud-e-b1-dft; then
  git switch cloud-e-b1-dft
else
  REMOTE_E_HEAD=$(git rev-parse refs/remotes/origin/cloud-e-b1-dft)
  git switch -c cloud-e-b1-dft "$REMOTE_E_HEAD"
fi
git pull --ff-only origin cloud-e-b1-dft

# Use the existing GPAW environment; do not install or upgrade packages here.
source /workspace/.venvs/gpaw-mpi/bin/activate
python -c "import gpaw, ase; print('GPAW', gpaw.__version__, 'ASE', ase.__version__)"

TASK_PACKAGE=research/mace-v12-transfer/coordination/reports/window-e/tsc_backbone_dft_batch_20261017
cd "$TASK_PACKAGE"
sha256sum -c SHA256SUMS.txt

# Non-scientific four-rank CLI smoke test; --help must exit 0 and creates no calculation output.
mpiexec -n 4 gpaw python -- "$PWD/run_surface_singlepoint.py" --help

# Preserve the original pre-SCF launcher failure at results/<TASK_ID>/.
# Retry it only at the distinct attempt-0002 path.
TASK_ID=TSC100_TI0001_COD9009647_2x2x2_REF
ATTEMPT_ID=attempt-0002
LOCK_SHA=579f5803cd6c3e5d522ecf396c4f25fb06ea9fe8939f9ef684e750702f9edb46
python execute_one_surface_task.py --task-id "$TASK_ID" --attempt-id "$ATTEMPT_ID" --lock-sha256 "$LOCK_SHA" --confirm-5h-window --preflight-only
# Run only if output status is PREFLIGHT_PASS.
python execute_one_surface_task.py --task-id "$TASK_ID" --attempt-id "$ATTEMPT_ID" --lock-sha256 "$LOCK_SHA" --confirm-5h-window

# Each completed task is archived and SHA-verified from a fresh fetch before proceeding.
# For each remaining task, omit --attempt-id and use the new LOCK_SHA above:
# TSC100_TI0001_COD9009647_2x2x2_RESP_Z020
# TSC100_SI0001_COD9009647_2x2x2_REF
# TSC100_SI0001_COD9009647_2x2x2_RESP_Z020
```

不要同时启动多个任务。每项固定使用 GPAW 26.7.0 / ASE 3.29.0、PBE/PAW、PW500 eV、FD 0.10 eV、标准 MP `(6,6,1)`、MSR1、默认 Davidson niter=2、SCF maxiter=160、MPI=4、OMP=1。固定坐标和晶胞；不做优化。输出的能量与应力只作诊断，获批用途是力标签候选。执行器没有 checkpoint；时间、内存或磁盘停止时保留日志并归档失败状态，不重启、不重试。
