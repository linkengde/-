#!/usr/bin/env bash
set -euo pipefail
REPO="$(git rev-parse --show-toplevel)"
TRANSFER="$REPO/research/mace-v12-transfer"
SYNC="$TRANSFER/coordination/sync_tasks.py"
V14="$TRANSFER/periodic_interface_v4/mace_periodic_v14_interface_energy"
python3 "$SYNC" claim window-b
PY="${MACE_PYTHON:-}"
for candidate in /workspace/.venvs/mace-v12/bin/python /workspace/.venvs/mace-v14/bin/python /workspace/.venvs/mace-v14-assist/bin/python; do
  if [[ -z "$PY" && -x "$candidate" ]] && "$candidate" -c 'import ase, mace, torch; assert ase.__version__ == "3.29.0" and mace.__version__ == "0.3.16" and torch.__version__.startswith("2.5.1+cpu")' >/dev/null 2>&1; then PY="$candidate"; fi
done
[[ -n "$PY" ]] || { echo 'MACE CPU environment not prepared; run run_v13_force_localization.sh first.' >&2; exit 2; }
[[ -f "$V14/data/dataset_manifest.json" ]] || { echo 'A has not published the verified v14 dataset yet; do not train.' >&2; exit 2; }
"$PY" - "$REPO" "$V14" <<'PY'
import hashlib, json, sys
from pathlib import Path
import numpy as np
from ase.io import read
repo, v14 = map(Path, sys.argv[1:])
project = repo / 'research/mace-v12-transfer/periodic_interface_v4'
m = json.loads((v14/'data/dataset_manifest.json').read_text())
assert m['split_sizes'] == {'train':31,'valid':2,'test':3}, m['split_sizes']
assert m['energy_composition_matrix']['rank'] == 4
assert m['frozen_v13_validation_and_test_byte_identical'] is True
for name, expected in m['output_sha256'].items():
    assert hashlib.sha256((v14/f'data/{name}.extxyz').read_bytes()).hexdigest() == expected
assert len(m['reserved_blind_holdouts']) == 3
required = [
 (project/'pbe_interface_v14_main_holdouts','AgC_registry_holdout_v14_01'),
 (project/'pbe_interface_v14_main_holdouts','AgSi_registry_holdout_v14_01'),
 (project/'pbe_interface_v14_parallel_acquisition','AgTi_registry_holdout_v14_01'),
]
train = read(v14/'data/train.extxyz', index=':')
for folder, label in required:
    result=folder/'calculations'/label
    verify=json.loads((result/'verification.json').read_text())
    assert verify['status']=='PASS' and verify['role']=='v14_holdout', (label, verify)
    output=next(result.glob('*_PW_PBE.extxyz'))
    assert hashlib.sha256(output.read_bytes()).hexdigest()==verify['sha256'][output.name]
    held=read(output)
    for used in train:
        assert not (len(held)==len(used) and held.get_chemical_symbols()==used.get_chemical_symbols()
                    and np.allclose(held.cell.array,used.cell.array,atol=1e-12,rtol=0)
                    and np.allclose(held.positions,used.positions,atol=1e-10,rtol=0)), label
print('v14 data/holdout gates PASS; 31/2/3 splits, rank 4/4, all three blind labels excluded.')
PY
if [[ -f "$V14/checkpoints/MACE_periodic_v14_interface_energy_run-45.model" ]]; then
  echo 'v14 checkpoint already exists; inspect its evaluation. This runner will not retrain over it.' >&2
  exit 2
fi
export XDG_CACHE_HOME=/tmp/mace-v14-cache MPLCONFIGDIR=/tmp/mpl-v14-cache FC_CACHEDIR=/tmp/fontconfig-v14-cache
mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$FC_CACHEDIR"
python3 "$SYNC" progress window-b --job v14_training --state running --iteration 0 --note 'Starting the one authorized MACE v14 CPU training run.'
python3 "$TRANSFER/coordination/reports/window-b/watch_mace_training.py" --parent-pid "$$" >/tmp/window-b-v14-progress.log 2>&1 &
WATCH_PID=$!
cleanup() { kill "$WATCH_PID" 2>/dev/null || true; }
trap cleanup EXIT
if ! MACE_RUN_TRAIN="$(dirname "$PY")/mace_run_train" bash "$V14/run_v14_training.sh"; then
  python3 "$SYNC" progress window-b --job v14_training --state failed --iteration 0 --note 'v14 training command exited with an error; preserve outputs and inspect train.stdout.' || true
  exit 1
fi
[[ -s "$V14/checkpoints/MACE_periodic_v14_interface_energy_run-45.model" ]] || { echo 'Training exited without the selected v14 checkpoint.' >&2; exit 2; }
kill "$WATCH_PID" 2>/dev/null || true
python3 "$SYNC" progress window-b --job v14_training --state completed --iteration 80 --note 'MACE v14 training exited successfully; selected checkpoint exists.'
python3 "$SYNC" progress window-b --job v14_evaluation --state running --iteration 0 --note 'Scoring frozen test data and three reserved v14 holdouts.'
if ! "$PY" "$V14/evaluate_v14_cycle.py" v14 > "$V14/results/evaluate_v14_cycle_v14.stdout"; then
  python3 "$SYNC" progress window-b --job v14_evaluation --state failed --iteration 0 --note 'v14 evaluation exited with an error; preserve model and inspect stdout.' || true
  exit 1
fi
"$PY" - "$V14/results/v14_independent_holdout_comparison.json" <<'PY'
import json, sys
d=json.load(open(sys.argv[1]))
print('v14 evaluation status:',d['overall_status'])
for kind, models in d['by_interface'].items():
    row=models['v14']
    print(kind, 'energy_MAE_meV_atom=',round(row['energy_MAE_meV_atom'],4),
          'force_RMSE_eV_A=',round(row['force_vector_RMSE_eV_A'],5),
          'max_separating_force_error_eV_A=',round(row['separating_force_max_abs_error_eV_A'],5),
          'screen_pass=',row['screen_pass'])
PY
python3 "$SYNC" progress window-b --job v14_evaluation --state completed --iteration 1 --note 'v14 frozen-test/independent-holdout evaluation written; inspect per-interface screen before downstream use.'
python3 "$SYNC" publish window-b
