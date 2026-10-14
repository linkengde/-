# DFT Resource Estimate — Minimal Interface Designs

## Basis and limitation

The current environment's previous gate is `RESOURCE_GATE_PASS=NO`: a 32-GiB cgroup had about 5.62 GiB current headroom, a historical peak only 73,728 B below its limit, and 15 `memory.max` events. A 394-atom interface initialization exceeded a 32-GiB environment. No successful matched peak-memory measurement exists for the 217-atom interface, so neither the current model nor either reduced design has a demonstrated RAM requirement.

The estimates below are geometric and algorithmic planning indicators only. They are not GPAW benchmarks, do not authorize SCF, and must not be used to declare a candidate safe in 32 GiB.

## Relative size indicators

| Model | Atoms | In-plane area vs 217-atom parent | Approximate cell volume vs parent | Reciprocal-space plan for similar in-plane sampling density | Status |
|---|---:|---:|---:|---|---|
| Archived 217-atom interface | 217 | 1.000 | 1.000 | Existing 2×2×1 proposal; convergence not demonstrated | Existing inputs remain unchanged |
| Recommended M77 shared-strain | 77 | About 0.347 | About 0.347 with same c length | 4×4×1 as first comparable-density mesh; 5×5×1 as denser check | Design only; resource not measured |
| M77 with substrate fixed | 77 | 7/19 = 0.368 | About 0.368 with same c length | 3×3×1 as first comparable-density mesh; 4×4×1 as denser check | Design alternative; resource not measured |
| Compact M52 | 52 | 4/19 = 0.211 | About 0.237 if c is ~42.66 Å for the six-layer Ag film | 4×4×1 then 5×5×1; use 5×5×1 as the denser check | Design only; resource not measured |

The smaller cells reduce plane-wave basis size at fixed cutoff because the cell volume is smaller, while the increased k-point meshes recover approximately the parent model's reciprocal-space sampling density. A rough coefficient-storage proxy proportional to cell volume × atom count × k-point count is about 0.49 of the 217-atom/2×2×1 proxy for shared-strain M77/4×4×1, about 0.52 for substrate-fixed M77/4×4×1, and about 0.35 for M52/5×5×1. This proxy omits FFT work arrays, eigensolver workspace, parallel distribution, PAW augmentation, allocator overhead, and checkpoints; it is not a peak-RAM estimate.

Shared-strain M77 is preferred as a balance: its area-7 cell has a larger lateral repeat, retains four Ag layers, and limits the maximum nominal in-plane strain to about 3.1%. M52 has fewer atoms but stronger finite-size/periodic-image effects and uses six Ag layers to recover vertical film thickness.

## Expected workload and stop conditions

The minimum proposed geometric scan is three fixed-geometry single points for M77 at `d_ref−10%`, `d_ref`, and `d_ref+10%`. If the optional M52 fallback is selected, it also requires three points; do not run both families without a separate scope decision. Each scan point needs a separate fixed geometry and output directory.

Maintain the original approved method when a later execution is separately released: GPAW 26.7.0, ASE 3.29.0, PBE/PAW, plane-wave cutoff 500 eV, Fermi–Dirac width 0.10 eV, convergence tolerances 1e-5, Pulay beta 0.05 / `nmaxold=8` / weight 100, four MPI ranks and one OpenMP thread. K-point convergence for the newly sized cells must be demonstrated; the proposed meshes above remain unrun.

Before release, use an external high-memory allocation with task-specific measured headroom, persistent scratch/checkpoint capacity, and monitored peak RSS/cgroup use. The former 64-GiB planning suggestion was not a benchmark or sufficiency guarantee. Stop on memory pressure, missing PAW data, hash mismatch, nonconvergence, or inadequate archive capacity; do not lower cutoff, loosen convergence, or reduce k-point density to fit.

The 32-GiB current workspace is not cleared for these calculations. Smaller atom counts make a pilot plausible, but until coordinate files exist and a preflight plus measured pilot establishes safe resource use, record `DFT_RESOURCE_STATUS=UNVERIFIED` and keep SCF stopped.
