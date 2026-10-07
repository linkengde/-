#!/usr/bin/env bash
# Draft only. It intentionally refuses to launch until A supplies external authorization.
set -euo pipefail

PACK_ROOT=""
WORKSPACE_ROOT=""
INPUT_ROOT=""
RUN_ROOT=""
ARCHIVE_ROOT=""
AUTHORIZATION=""
OWNER=""
LABELS=()

while (($#)); do
  case "$1" in
    --pack-root) PACK_ROOT="${2:?missing value for --pack-root}"; shift 2 ;;
    --workspace-root) WORKSPACE_ROOT="${2:?missing value for --workspace-root}"; shift 2 ;;
    --input-root) INPUT_ROOT="${2:?missing value for --input-root}"; shift 2 ;;
    --run-root) RUN_ROOT="${2:?missing value for --run-root}"; shift 2 ;;
    --archive-root) ARCHIVE_ROOT="${2:?missing value for --archive-root}"; shift 2 ;;
    --authorization) AUTHORIZATION="${2:?missing value for --authorization}"; shift 2 ;;
    --owner) OWNER="${2:?missing value for --owner}"; shift 2 ;;
    --labels)
      shift
      while (($#)) && [[ "$1" != --* ]]; do LABELS+=("$1"); shift; done
      ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

for pair in \
  "pack-root:$PACK_ROOT" "workspace-root:$WORKSPACE_ROOT" "input-root:$INPUT_ROOT" \
  "run-root:$RUN_ROOT" "archive-root:$ARCHIVE_ROOT" "authorization:$AUTHORIZATION" "owner:$OWNER"; do
  [[ "${pair#*:}" != "" ]] || { echo "Missing required --${pair%%:*}" >&2; exit 2; }
done
((${#LABELS[@]} > 0)) || { echo "Supply one or more --labels" >&2; exit 2; }
[[ "$OWNER" == window-a || "$OWNER" == window-b ]] || { echo "Invalid owner" >&2; exit 2; }

PACK_ROOT="$(realpath "$PACK_ROOT")"
WORKSPACE_ROOT="$(realpath "$WORKSPACE_ROOT")"
INPUT_ROOT="$(realpath "$INPUT_ROOT")"
RUN_ROOT="$(realpath -m "$RUN_ROOT")"
ARCHIVE_ROOT="$(realpath -m "$ARCHIVE_ROOT")"
AUTHORIZATION="$(realpath "$AUTHORIZATION")"

GPAW_PYTHON="${GPAW_PYTHON:-/workspace/.venvs/gpaw-mpi/bin/python}"
[[ -x "$GPAW_PYTHON" ]] || { echo "GPAW Python is unavailable: $GPAW_PYTHON" >&2; exit 2; }
GPAW_CLI="$(dirname "$GPAW_PYTHON")/gpaw"
[[ -x "$GPAW_CLI" ]] || { echo "GPAW CLI is unavailable: $GPAW_CLI" >&2; exit 2; }

# This must pass before GPAW import or MPI calculation; use the pinned env for ASE.
"$GPAW_PYTHON" "$PACK_ROOT/execution_preflight.py" \
  --workspace-root "$WORKSPACE_ROOT" --pack-root "$PACK_ROOT" --input-root "$INPUT_ROOT" \
  --run-root "$RUN_ROOT" --archive-root "$ARCHIVE_ROOT" --authorization "$AUTHORIZATION" \
  --owner "$OWNER" --labels "${LABELS[@]}"

mkdir -p "$RUN_ROOT" "$ARCHIVE_ROOT"
exec 9>"$RUN_ROOT/.v17_holdout_execution.lock"
flock -n 9 || { echo "Another V17 holdout runner holds the execution lock" >&2; exit 2; }
if pgrep -f '[d]ft_single_point.py' >/dev/null || pgrep -x mpirun >/dev/null || pgrep -x prterun >/dev/null; then
  echo "An MPI/DFT job is already active in this environment; refusing overlap" >&2
  exit 2
fi

MPI_DEPS="${MPI_DEPS:-/workspace/.local/gpaw-mpi-deps}"
export GPAW_MPI_BACKEND=cgpaw
export LD_LIBRARY_PATH="$MPI_DEPS/usr/lib/x86_64-linux-gnu/openmpi/lib:$MPI_DEPS/usr/lib/x86_64-linux-gnu:/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu:/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export GPAW_SETUP_PATH="$("$GPAW_PYTHON" -c 'import gpaw_data; print(gpaw_data.datapath())')"

"$GPAW_PYTHON" - <<'PY'
import os
from pathlib import Path
available = len(os.sched_getaffinity(0))
quota_file = Path('/sys/fs/cgroup/cpu.max')
if quota_file.exists():
    quota, period = quota_file.read_text().split()
    if quota != 'max':
        available = min(available, int(quota) // int(period))
if available < 4:
    raise SystemExit('Four actual CPU cores are required; refusing to oversubscribe.')
import gpaw
if gpaw.__version__ != '26.7.0':
    raise SystemExit(f'Expected GPAW 26.7.0; found {gpaw.__version__}')
PY
"$(dirname "$GPAW_PYTHON")/gpaw" info | rg -q 'MPI enabled +yes' || { echo "GPAW MPI support missing" >&2; exit 2; }
"$(dirname "$GPAW_PYTHON")/gpaw" info | rg -q 'scalapack +yes' || { echo "GPAW ScaLAPACK support missing" >&2; exit 2; }

for LABEL in "${LABELS[@]}"; do
  # Repeat the non-overwrite and authorization checks immediately before each label.
  "$GPAW_PYTHON" "$PACK_ROOT/execution_preflight.py" \
    --workspace-root "$WORKSPACE_ROOT" --pack-root "$PACK_ROOT" --input-root "$INPUT_ROOT" \
    --run-root "$RUN_ROOT" --archive-root "$ARCHIVE_ROOT" --authorization "$AUTHORIZATION" \
    --owner "$OWNER" --labels "$LABEL"
  mpirun --bind-to core --map-by core -n 4 "$GPAW_CLI" python "$PACK_ROOT/dft_single_point.py" \
    --workspace-root "$WORKSPACE_ROOT" --pack-root "$PACK_ROOT" --input-root "$INPUT_ROOT" \
    --run-root "$RUN_ROOT" --authorization "$AUTHORIZATION" --owner "$OWNER" --label "$LABEL" \
    >"$RUN_ROOT/$LABEL.launcher.log" 2>&1
  python3 "$PACK_ROOT/archive_verify.py" \
    --workspace-root "$WORKSPACE_ROOT" --pack-root "$PACK_ROOT" --input-root "$INPUT_ROOT" \
    --run-root "$RUN_ROOT" --archive-root "$ARCHIVE_ROOT" --authorization "$AUTHORIZATION" \
    --owner "$OWNER" --label "$LABEL"
done

echo "Assigned V17 holdout labels completed and compact archives verified."
