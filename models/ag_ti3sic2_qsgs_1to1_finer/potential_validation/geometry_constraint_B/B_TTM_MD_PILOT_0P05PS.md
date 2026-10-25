# Model B: reference TTM–MD coupling pilot (0.05 ps)

Date: 2026-10-25. This is a bounded software/physics pilot on the geometry-constrained exploratory Model B. It is not a production TTM implementation, not a cold-baseline acceptance, and not a melt-pool result.

## Setup

- Geometry: the 916,534-atom relaxed-parent QSGS candidate; the prior Ag-only FIRE-100 coordinates are loaded.
- Potential: Ag EAM, Ti–Si–C TSCBOP, and the archived WCA-only geometric obstacle. Ti/Si/C and the bottom layer are fixed; Ag above 20 Å is mobile in NVE.
- Source/coupling: the repository Python reference TTM solver, 10 Å nominal electron grid, `R=30 Å`, `l_skin=20 Å`, `I0=5 eV/(ps·Å²)`, and the candidate Ag parameter table. The opt-in sub-300 K continuation was enabled (`Ce∝T`, constant `G` and `ke` below 300 K). The source and electron–lattice exchange act on Newton-layer Ag. The Ag-cell map is held at its initial geometry during this short pilot.
- Integration: 0.5 fs MD step, 100 LAMMPS MD steps (0.05 ps). The Python driver advances LAMMPS by one step, updates the TTM grid, and writes the equal/opposite lattice energy transfer back to Ag velocities. Each MD step required 242–351 explicit electron substeps.
- No liquid-cluster classifier or automatic heat-off controller was connected in this pilot.

## Results

| Time | Newton Ag kinetic temperature | Max active-cell electron temperature | Cumulative source energy | Cumulative electron-to-lattice transfer | Combined energy residual* |
|---:|---:|---:|---:|---:|---:|
| 0.0005 ps | 299.58 K | 581.86 K | 1.467 eV | 0.000496 eV | −0.000005 eV |
| 0.0100 ps | 291.45 K | 822.91 K | 29.347 eV | 0.172 eV | 0.159 eV |
| 0.0200 ps | 267.94 K | 871.57 K | 58.695 eV | 0.720 eV | 0.701 eV |
| 0.0300 ps | 233.17 K | 913.23 K | 88.042 eV | 1.650 eV | 1.636 eV |
| 0.0400 ps | 193.99 K | 950.94 K | 117.390 eV | 2.968 eV | 2.966 eV |
| 0.0500 ps | 158.95 K | 985.64 K | 146.737 eV | 4.644 eV | 4.653 eV |

\* Residual is total LAMMPS energy plus electron-grid energy and transferred lattice energy, minus the initial combined energy and cumulative source input. The final `run 0` recalculation gave 4.745 eV residual, about 3.24% of injected energy. The TTM-only one-step ledger closed to `7.3×10⁻¹⁵ eV`; the larger combined residual accumulates during MD and should be reduced before any longer trajectory is interpreted quantitatively.

The short run completed without a TTM table-range exception. It did **not** form a stable liquid Ag cluster: the Ag kinetic temperature fell from 300 K to 159 K while electronic temperature rose. Only 4.64 eV of 146.74 eV injected energy reached the Newton-layer Ag lattice over 0.05 ps. This is consistent with continued structural/energy relaxation dominating the short pilot; it is not evidence that B cannot melt under a longer, validated run.

## Decision

The API bridge now demonstrates an actual LAMMPS MD step followed by the reference TTM update and Ag velocity energy handoff. That is a useful integration milestone. It does not clear either outstanding gate:

1. The B model still lacks an accepted cold starting state; its no-source 100 fs screen previously ended at 183.98 K and showed 63 Ag atoms outside the largest connectivity component.
2. The Python reference is far too costly for production-length simulation at this grid and timestep, and the accumulated MD-plus-TTM energy residual needs investigation. A compiled/accelerated implementation and longer real-MD energy validation are still required before the 100 ps target.

No stable-pool shutoff, 100 ps run, or cooling/solidification run was performed. Do not interpret this as a physical prediction of Ag–Ti₃SiC₂ melting or interface heat transfer; the fixed scaffold and repulsive WCA cross interaction are geometry-only assumptions.

## Reproduction

The LAMMPS Python library in this environment links to MPICH (`libmpi.so.12`), while the system `mpi4py` is linked to OpenMPI. The recorded run used a temporary `mpi4py` build against the repository MPICH compiler wrapper:

```bash
UV_CACHE_DIR=/tmp/ag-ttm-uv-cache \
MPICC=/workspace/.venvs/lammps/bin/mpicc \
uv pip install --target /tmp/ag-ttm-mpich --no-binary mpi4py mpi4py

LD_LIBRARY_PATH=/workspace/.venvs/lammps/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} \
PYTHONPATH=/tmp/ag-ttm-mpich:/workspace/.venvs/lammps/lib/python3.12/site-packages \
OMP_NUM_THREADS=4 python \
models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/geometry_constraint_B/b_ttm_reference_md_pilot.py
```

Run from the repository root. The driver and captured output are `b_ttm_reference_md_pilot.py` and `log.b_ttm_reference_md_pilot`; `SHA256SUMS_TTM_PILOT.txt` records their hashes.
