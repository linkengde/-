#!/workspace/venvs/gpaw/bin/python
"""Run one frozen V13 validation geometry through the training DFT protocol."""
import argparse
import importlib.util
import json
from pathlib import Path

WORK = Path('/workspace/potential-fit/round7_pod_snap/models/target_200x200x300_qsgs_material200/v13_source_and_fit_20261123')
CANDIDATES = WORK.parent / 'v13_candidate_sources_20261123'
MANIFEST = CANDIDATES / 'validation_candidate_manifest_v13.json'
LABEL_ROOT = CANDIDATES / 'validation_dft_labels'
TRAIN_RUNNER = WORK / 'run_v13_candidate_dft.py'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--candidate-id', required=True)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('v13_train_dft_runner', TRAIN_RUNNER)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    runner.LABEL_ROOT = LABEL_ROOT
    manifest = json.loads(MANIFEST.read_text())
    rows = [r for r in manifest['candidates'] if r['candidate_id'] == args.candidate_id]
    if len(rows) != 1:
        raise SystemExit(f"validation candidate {args.candidate_id} not found exactly once")
    result = runner.run_one(rows[0])
    if result.get('run_status') != 'complete':
        raise SystemExit(f"validation DFT did not produce a complete label: {result.get('run_status')}")


if __name__ == '__main__':
    main()
