#!/usr/bin/env bash
set -euo pipefail
REPO="$(git rev-parse --show-toplevel)"
REPORT="$REPO/research/mace-v12-transfer/coordination/reports/window-b"
SYNC="$REPO/research/mace-v12-transfer/coordination/sync_tasks.py"
python3 "$SYNC" claim window-b
PY="${MACE_PYTHON:-}"
for candidate in /workspace/.venvs/mace-v12/bin/python /workspace/.venvs/mace-v14/bin/python /workspace/.venvs/mace-v14-assist/bin/python; do
  if [[ -z "$PY" && -x "$candidate" ]] && "$candidate" -c 'import ase, mace, torch; assert ase.__version__ == "3.29.0" and mace.__version__ == "0.3.16" and torch.__version__.startswith("2.5.1+cpu")' >/dev/null 2>&1; then PY="$candidate"; fi
done
if [[ -z "$PY" ]]; then
  command -v python3.12 >/dev/null || { echo 'Python 3.12 is required; see environment_setup.md.' >&2; exit 2; }
  PY=/workspace/.venvs/mace-v14-assist/bin/python
  if [[ ! -x "$PY" ]]; then
    python3.12 -m venv /workspace/.venvs/mace-v14-assist
    PIP_CACHE_DIR=/tmp/pip-cache /workspace/.venvs/mace-v14-assist/bin/python -m pip install --upgrade pip
    PIP_CACHE_DIR=/tmp/pip-cache /workspace/.venvs/mace-v14-assist/bin/python -m pip install --index-url https://download.pytorch.org/whl/cpu 'torch==2.5.1+cpu'
    PIP_CACHE_DIR=/tmp/pip-cache /workspace/.venvs/mace-v14-assist/bin/python -m pip install 'mace-torch==0.3.16' 'ase==3.29.0'
  fi
  PY=/workspace/.venvs/mace-v14-assist/bin/python
fi
"$PY" -c 'import ase, mace, torch; print("ASE",ase.__version__,"MACE",mace.__version__,"PyTorch",torch.__version__)'
export XDG_CACHE_HOME=/tmp/mace-v14-cache MPLCONFIGDIR=/tmp/mpl-v14-cache FC_CACHEDIR=/tmp/fontconfig-v14-cache
mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$FC_CACHEDIR"
"$(dirname "$PY")/mace_run_train" --help >/tmp/window-b-mace-run-train-help.txt
python3 "$SYNC" progress window-b --job v13_force_localization --state running --iteration 0 --note 'MACE environment ready; localizing v13 per-atom force residuals.'
"$PY" "$REPORT/analyze_v13_force_localization.py"
python3 "$SYNC" progress window-b --job v13_force_localization --state completed --iteration 1 --note 'Per-atom force residual report written; v14 holdouts and training data unchanged.'
python3 "$SYNC" publish window-b
