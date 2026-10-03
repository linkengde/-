# Environment setup for the v12 handoff

The tested original cloud uses Python 3.12 and two CPU virtual environments. The GPAW PAW setups are installed through `gpaw-data`; the DFT runner also sets `GPAW_SETUP_PATH` from `gpaw_data.datapath()`.

## GPAW single-point calculations

```bash
python3.12 -m venv .venv-gpaw
. .venv-gpaw/bin/activate
PIP_CACHE_DIR=/tmp/pip-cache python -m pip install \
  'gpaw==26.7.0' 'gpaw-data==1.2.1' 'ase==3.29.0'
python -c 'import ase, gpaw, gpaw_data; print(gpaw.__version__, ase.__version__, gpaw_data.datapath())'
```

To run the queued AgSi 2.2 A point from the package:

```bash
cd research/mace-v12-transfer/periodic_interface_v4/pbe_interface_energy_additions_v12
GPAW_PYTHON=/path/to/.venv-gpaw/bin/python ./run_AgSi_d2p2.sh
```

The run writes rolling state under `/tmp/pbe_interface_v12_runs/AgSi_d2p2` and archives compact converged outputs under `calculations/AgSi_d2p2/`. Do not copy a live `.gpw` file while GPAW is updating it.

## MACE inference and training

```bash
python3.12 -m venv .venv-mace
. .venv-mace/bin/activate
PIP_CACHE_DIR=/tmp/pip-cache python -m pip install \
  --index-url https://download.pytorch.org/whl/cpu 'torch==2.5.1+cpu'
PIP_CACHE_DIR=/tmp/pip-cache python -m pip install \
  'mace-torch==0.3.16' 'ase==3.29.0'
mace_run_train --help
```

The MACE-MP-0b3-medium foundation model is bundled under `periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/`; verify it against `file_manifest.json`. After all four DFT result archives are present and pass `archive_completed_pw.py`, build the dataset and run `./run_v12_training.sh` from the v12 directory.

The tested versions and current validation results are in `environment_validation.md`. The training command check does not train a model; v12 fitting remains blocked until all four PW-PBE references converge and are archived.
