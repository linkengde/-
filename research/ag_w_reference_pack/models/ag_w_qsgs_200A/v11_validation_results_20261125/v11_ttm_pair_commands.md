# Reproducing the V11 TTM software-link pair

Run in the repository root with the Python LAMMPS binding, SciPy, and NumPy installed. Set `LD_LIBRARY_PATH` so the Python binding can locate the MPI library paired with that LAMMPS build. The paths below use this checkout's V11 branch.

```bash
ROOT="$PWD"
BASE="$ROOT/research/ag_w_reference_pack/models/ag_w_qsgs_200A"
PILOT="$BASE/ttm/agw_ttm_md_pilot.py"
DATA="$BASE/v11_generated_data/agw_v11_gyroid.data.gz"
FRACTIONS="$BASE/v11_generated_data/agw_v11_geometry_fractions.npz"
MANIFEST="$BASE/v11_validation_results_20261125/v11_ttm_exploratory_manifest.json"
RESULTS="$BASE/v11_validation_results_20261125"
TMPDIR_TTM="$(mktemp -d)"
gzip -dc "$RESULTS/v11_cold.dump.gz" > "$TMPDIR_TTM/v11_cold.dump"
```

Use the same decompressed step-8200 restart for each run. The first command repeats the no-source control; the second applies the Ag Gaussian source. These are short software-link checks and use an explicit exploratory override.

```bash
python "$PILOT" \
  --data "$DATA" --audit "$MANIFEST" \
  --restart-dump "$TMPDIR_TTM/v11_cold.dump" --restart-step 8200 \
  --geometry-fractions-npz "$FRACTIONS" \
  --output "$TMPDIR_TTM/no_source.json" \
  --pulse-fwhm-ps 0.01 --pulse-center-ps 0.02 --cooling-ps 0.02 \
  --grid-A 10 --box-xy-A 200 --dt-ps 0.0001 \
  --source-scale 0 --initial-electron-temperature-from-lattice --exploratory

python "$PILOT" \
  --data "$DATA" --audit "$MANIFEST" \
  --restart-dump "$TMPDIR_TTM/v11_cold.dump" --restart-step 8200 \
  --geometry-fractions-npz "$FRACTIONS" \
  --output "$TMPDIR_TTM/gaussian_source.json" \
  --pulse-fwhm-ps 0.01 --pulse-center-ps 0.02 --cooling-ps 0.02 \
  --grid-A 10 --box-xy-A 200 --dt-ps 0.0001 \
  --source-scale 1 --initial-electron-temperature-from-lattice --exploratory
```

The first local run used the MPICH-compatible LAMMPS Python library by setting `LD_LIBRARY_PATH=/workspace/.venvs/ag-qsgs/lib:/workspace/.venvs/ag-qsgs/lib/python3.12/site-packages/lammps`; adjust this for another installation.
