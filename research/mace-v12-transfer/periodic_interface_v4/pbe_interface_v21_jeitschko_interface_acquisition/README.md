# Jeitschko–Nowotny 1967 Ag/Ti₃SiC₂ interface acquisition

This task builds three fixed-geometry, unlabeled 52-atom Ag/Ti₃SiC₂ slab inputs from the reconstructed 1967 single-crystal source. It is intended to acquire missing Ag–C, Ag–Si and Ag–Ti interface forces after the V19 development gate failed. It does not change V19 data, model weights, split assignments or training status.

## Source and geometry

The substrate derives from `Ti3SiC2_Jeitschko1967_fulltext_Fourier_variant.cif` (DOI `10.1007/BF00899949`), a full-text reconstruction rather than an ICSD 25762 export. The source CIF and uploaded paper hashes are in `input_manifest.json`. The model assumes full occupancy despite the paper's discussion of possible nonstoichiometry. Each input uses a 2×2 in-plane substrate repeat and one Ag monolayer with fixed fractional registry `(0.37, 0.21)`. The Ag plane sits 2.4 Å above the C-terminated slab or 2.6 Å above the Si/Ti-terminated slab. These are acquisition geometries, not measured or relaxed equilibrium interfaces.

`build_inputs.py` reads no energy or force fields. It checks the three generated files, screens geometry-only columns against the frozen V19 train36/dev3, residual proposal inputs and V17 development geometries, and records exact file hashes. The direct same-composition geometry screen found no exact or near duplicates. This is not a universal supercell-equivalence audit of every historical archive. The three interfaces share one 1967 parent family, so they do not constitute three independent source families.

## Method and role

Candidate calculations use GPAW 26.7.0 / ASE 3.29.0, PBE, 500 eV plane waves, Fermi–Dirac width 0.10 eV, 6×6×1 k points, four MPI ranks and an xy dipole correction. Native extrapolated energy, finite-width free energy and forces are stored separately. Training integration remains disabled until the energy target and source-family audit are decided.

The first registered run is the AgSi 6×6×1 candidate. A same-geometry 8×8×1 job is control-only. The provisional mesh budget is ±2 meV/atom for each energy, force-vector RMS ≤0.01 eV/Å and marked Ag–Si relative-force projection ≤0.02 eV/Å. The 30 GiB cgroup stop guard leaves 2 GiB below the 32 GiB hard limit; jobs run sequentially and are never automatically retried. `.gpw` state files are removed from archived results.

## Execution

Run from the repository root in the GPAW MPI environment:

```text
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_interface_acquisition/build_inputs.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_interface_acquisition/launch_preflight.py AgSi_JN1967_registryB_k6x6x1_sigma0p10_preflight_v21
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_interface_acquisition/launch_preflight.py AgSi_JN1967_registryB_k8x8x1_sigma0p10_control_v21
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_interface_acquisition/compare_kmesh.py
```

`launch_preflight.py` validates cgroup and disk headroom, MPI/ScaLAPACK, package and PAW hashes, and the shared DFT lock before each run. The runner preserves per-iteration progress, checks memory every 10 seconds, and stores result text plus extxyz and SHA256 inventory. A failed run is retained under `calculations/failed_attempts/` for diagnosis; do not overwrite it or retry without correcting the cause.

## Current status

Geometry candidates and method scripts are prepared. The first attempt stopped before SCF because the runner passed a `PoissonSolver` object where GPAW 26.7 expects the `{"dipolelayer": "xy"}` parameter mapping. No energy or force label was produced. That attempt is preserved with its launcher log and classified as a runner initialization error; the runner was corrected before the next attempt. See `calculations/` for live or completed status and `input_manifest.json` for job flags.

Even if the AgSi mesh comparison passes, only the same 1967-source AgC, AgSi and AgTi geometries should proceed as acquisition candidates; the 8×8×1 endpoint remains a numerical control. A successful calculation does not by itself resolve V19's failed development gate or grant a train/dev role.
