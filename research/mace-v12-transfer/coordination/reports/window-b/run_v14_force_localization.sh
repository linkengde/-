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
[[ -n "$PY" ]] || { echo 'Prepared MACE CPU environment not found; do not reinstall or retrain v14. Check the B environment.' >&2; exit 2; }
"$PY" -c 'import ase,mace,torch; print("ASE",ase.__version__,"MACE",mace.__version__,"PyTorch",torch.__version__)'
export XDG_CACHE_HOME=/tmp/mace-v14-cache MPLCONFIGDIR=/tmp/mpl-v14-cache FC_CACHEDIR=/tmp/fontconfig-v14-cache
mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$FC_CACHEDIR"
python3 "$SYNC" progress window-b --job v14_force_localization --state running --iteration 0 --note 'Post-hoc localization of the three published v14 holdout force failures.'
"$PY" "$REPORT/analyze_v14_force_localization.py"
python3 "$SYNC" progress window-b --job v14_force_localization --state completed --iteration 1 --note 'Per-atom v14 holdout residual report written; no training or blind data changed.'
python3 "$SYNC" publish window-b
