# Ag/Ti₃SiC₂ transient FVM solver optimization

This folder adds an optional Numba CPU path and instruments the existing single-temperature enthalpy solver. The reference Python/SciPy solver and `phase_change_validation/` inputs are left intact. No arc boundary, liquid flow, Marangoni stress, arc pressure, evaporation, or MD changes are included.

## Environment and install

The measured cloud allocation was 4 CPU cores, 32 GiB memory, no GPU. NumPy 2.3.5 and SciPy 1.17.0 were available. Numba 0.67.0 was initially absent and was installed under `/tmp/solver-opt-deps` for this session; it was not installed globally. For a repeatable environment, install the optional dependency:

```bash
python -m pip install -r models/ag_ti3sic2_qsgs_ultrafine_micro_fvm/solver_optimization/requirements-numba.txt
```

The Numba path currently supports uniform cubic 3-D grids. The pure Python path remains the reference and does not require Numba. PyAMG is an optional pilot dependency (`requirements-multigrid.txt`); its small hot-state probe was slower than FFT after excluding JIT warm-up, so the production benchmark defaults to FFT.

## Reproduce

From the repository root:

```bash
export PYTHONPATH="models/ag_ti3sic2_qsgs_ultrafine_micro_fvm/phase_change_validation:models/ag_ti3sic2_qsgs_ultrafine_micro_fvm/solver_optimization${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m unittest discover -s models/ag_ti3sic2_qsgs_ultrafine_micro_fvm/solver_optimization -p 'test_*.py' -v
python models/ag_ti3sic2_qsgs_ultrafine_micro_fvm/solver_optimization/summarize_results.py
python models/ag_ti3sic2_qsgs_ultrafine_micro_fvm/solver_optimization/plot_results.py
```

Example benchmark commands (use a new output directory for a new run):

```bash
python models/ag_ti3sic2_qsgs_ultrafine_micro_fvm/solver_optimization/run_benchmark.py \
  --case bicontinuous_qsgs_seed1105 --n 50 --solver numba \
  --linear-policy ew_safe --preconditioner fft --threads 4 --dt-factor 0.25 \
  --checkpoint-dir models/ag_ti3sic2_qsgs_ultrafine_micro_fvm/solver_optimization/results/checkpoints/my_run
```

The benchmark writes a flushed per-step JSONL log and an atomic float64 temperature/enthalpy checkpoint after every accepted step. `--resume` resumes that run only if the geometry hash, source hash, material settings, solver tolerances, and run configuration match. Static Ag fraction and source-rate fields are shared by resolution under `results/checkpoints/shared_inputs/` and referenced from each run with relative symlinks. `--max-steps` is useful for short probes; `--wall-limit-s` stops at a recoverable accepted-step boundary. Thread counts are explicitly limited to 1, 2, or 4; BLAS/OMP/MKL thread counts are fixed at 1 to avoid oversubscription and remap nondeterminism. The `ew_safe` inexact-Newton policy caps early PCG forcing at 0.1 and keeps the nonlinear acceptance tolerance unchanged; `ew` is retained as a less conservative comparison.

The optional AMG comparison is a single preheated 50³ step, not a cycle:

```bash
python -m pip install -r models/ag_ti3sic2_qsgs_ultrafine_micro_fvm/solver_optimization/requirements-multigrid.txt
python models/ag_ti3sic2_qsgs_ultrafine_micro_fvm/solver_optimization/run_amg_probe.py
```

## Outputs

- `SOLVER_PROFILING.md`: per-step/iteration bottlenecks and CPU profile.
- `GEOMETRY_MAPPING.md`: conservative mapping, mixed cells, source-integral error, and connectivity caveats.
- `OPTIMIZATION_COMPARISON.csv`: machine-readable run comparisons.
- `NONLINEAR_CONVERGENCE.md`: PCG forcing, Newton, retry, and AMG results.
- `CHECKPOINT_VALIDATION.md`: persistence and bitwise restart check.
- `NUMERICAL_REGRESSION.md`: reference/optimized kernel and phase-change checks.
- `MATERIAL_PROPERTY_AUDIT.md`: prior candidate-property audit and unchanged NO-GO items.
- `PERFORMANCE_REPORT.md` and `FVM_READINESS.md`: measured stage decisions.
- `results/checkpoints/`: selected full and partial restart states with per-step traces.
- `figures/`: timing and threshold-sensitivity plots.

To verify all saved scripts, reports, CSVs, figures, logs, and checkpoint files, run from this folder:

```bash
python write_manifest.py
sha256sum -c SHA256SUMS.txt
```

The volumetric source is a synthetic localized numerical pulse. It is not a measured arc heat boundary. All phase-change numbers use the existing candidate properties and remain exploratory until the material-source audit is resolved.
