# Jeitschko–Nowotny 1967 Ag/Ti₃SiC₂ interface acquisition

This task builds three fixed-geometry, unlabeled 52-atom Ag/Ti₃SiC₂ slab inputs from the reconstructed 1967 single-crystal source. It is intended to acquire missing Ag–C, Ag–Si and Ag–Ti interface forces after the V19 development gate failed. It does not change V19 data, model weights, split assignments or training status.

## Source and geometry

The substrate derives from `Ti3SiC2_Jeitschko1967_fulltext_Fourier_variant.cif` (DOI `10.1007/BF00899949`), a full-text reconstruction rather than an ICSD 25762 export. The source CIF and uploaded paper hashes are in `input_manifest.json`. The model assumes full occupancy despite the paper's discussion of possible nonstoichiometry. Each input uses a 2×2 in-plane substrate repeat and one Ag monolayer with fixed fractional registry `(0.37, 0.21)`. The Ag plane sits 2.4 Å above the C-terminated slab or 2.6 Å above the Si/Ti-terminated slab. These are acquisition geometries, not measured or relaxed equilibrium interfaces.

`build_inputs.py` reads no energy or force fields. It checks the three generated files, screens geometry-only columns against the frozen V19 train36/dev3, residual proposal inputs and V17 development geometries, and records exact file hashes. The direct same-composition geometry screen found no exact or near duplicates. This is not a universal supercell-equivalence audit of every historical archive. The three interfaces share one 1967 parent family, so they do not constitute three independent source families.

## Method and role

Candidate calculations use GPAW 26.7.0 / ASE 3.29.0, PBE, 500 eV plane waves, Fermi–Dirac width 0.10 eV, 6×6×1 k points, four MPI ranks and an xy dipole correction. Native extrapolated energy, finite-width free energy and forces are stored separately. Training integration remains disabled until the energy target and source-family audit are decided.

The first registered run is the AgSi 6×6×1 candidate. A same-geometry 8×8×1 job is control-only. The provisional mesh budget is ±2 meV/atom for each energy, force-vector RMS ≤0.01 eV/Å and marked Ag–Si relative-force projection ≤0.02 eV/Å. The 30 GiB cgroup stop guard leaves 2 GiB below the 32 GiB hard limit; jobs run sequentially and are never automatically retried. `.gpw` state files are removed from archived results.

## Execution

Run from the repository root in the GPAW MPI environment. `run_mesh_sequence.py` waits for the already-running AgSi 6×6×1 job, verifies it, runs the same-geometry 8×8×1 control and compares both meshes. It starts the AgC and AgTi 6×6×1 labels sequentially only if the comparison passes its provisional budgets:

```text
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_interface_acquisition/build_inputs.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_interface_acquisition/run_mesh_sequence.py
```

`launch_preflight.py` validates cgroup and disk headroom, MPI/ScaLAPACK, package and PAW hashes, and the shared DFT lock before each run. The runner preserves per-iteration progress, checks memory every 10 seconds, and stores result text plus extxyz and SHA256 inventory. A failed run is retained under `calculations/failed_attempts/` for diagnosis; do not overwrite it or retry without correcting the cause.

## Current status

The corrected AgSi 6×6×1 fixed-geometry calculation converged in 59 SCF iterations. It produced native extrapolated energy `-390.0466850492484 eV/cell`, finite-width free energy `-390.44621477785967 eV/cell`, maximum force `1.55115914 eV/Å` and force-vector RMS `0.77551266 eV/Å`. These are single-point labels on the registered unrelaxed geometry; the nonzero forces do not indicate SCF failure. The result archive is `calculations/AgSi_JN1967_registryB_k6x6x1_sigma0p10_preflight_v21/` and passes its integrity and method checks. Its recovery note records that the scientific runner wrote complete outputs but the outer launcher did not write its runtime record or finish packaging; the unavailable peak-memory value is left unknown.

The first attempt stopped before SCF because the runner passed a `PoissonSolver` object where GPAW 26.7 expects the `{"dipolelayer": "xy"}` parameter mapping. No energy or force label was produced. That attempt is preserved with its launcher log and classified as a runner initialization error; the runner was corrected before the successful run.

The same-geometry AgSi 8×8×1 control converged in 56 SCF iterations. Relative to 6×6×1, the native-energy difference is `-0.1705 meV/atom`, free-energy difference is `-0.1213 meV/atom`, force-vector RMS difference is `0.00298 eV/Å`, and marked Ag–Si pair projection difference is `0.0000882 eV/Å`. All are within the provisional mesh budgets. The control remains numerical-only and is not a training endpoint. Its archived peak memory was 13.96 GiB. The full result and comparison are in `calculations/AgSi_JN1967_registryB_k8x8x1_sigma0p10_control_v21/` and `comparisons/comparison_AgSi_JN1967_registryB_k6x6x1_sigma0p10_preflight_v21__AgSi_JN1967_registryB_k8x8x1_sigma0p10_control_v21.json`.

An earlier 8×8 execution was interrupted when its workspace disappeared; no checkpoint or result from that execution was recovered. The archived control is the clean rerun from the registered input. A first archive verification also exposed an over-strict `1e-12 eV/Å` force-summary comparison against force components rounded in extxyz. The verifier now allows `1e-8 eV/Å`, accounting only for serialization roundoff.

The AgC 6×6×1 fixed-geometry candidate converged in 72 SCF iterations and passed archive verification. It produced native extrapolated energy `-378.6672443061853 eV/cell`, free energy `-379.1163275989818 eV/cell`, maximum force `2.97877290 eV/Å` and force-vector RMS `1.14779784 eV/Å`. These are single-point labels on the registered unrelaxed geometry. Peak memory was 9.94 GiB. Its role remains `train_acquisition_candidate_pending_audit`, with training integration disabled; see `calculations/AgC_JN1967_registryB_k6x6x1_sigma0p10_acquisition_v21/`.

The AgTi 6×6×1 fixed-geometry candidate converged in 49 SCF iterations and passed archive verification. It produced native extrapolated energy `-375.7532401992928 eV/cell`, free energy `-376.298068189782 eV/cell`, maximum force `1.80222951 eV/Å` and force-vector RMS `0.73076157 eV/Å`. These are single-point labels on the registered unrelaxed geometry. Peak observed memory was 9.27 GiB. Its role remains `train_acquisition_candidate_pending_audit`, with training integration disabled; the archive, launcher record and runtime record are in `calculations/AgTi_JN1967_registryB_k6x6x1_sigma0p10_acquisition_v21/`.

Only the same 1967-source AgC, AgSi and AgTi geometries proceed as acquisition candidates; the 8×8×1 endpoint remains a numerical control. The mesh comparison does not resolve V19's failed development gate or grant a train/dev role.
