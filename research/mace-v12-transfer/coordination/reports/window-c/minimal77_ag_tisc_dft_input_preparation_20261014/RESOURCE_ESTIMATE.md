# Resource estimate — 77-atom versus archived 217-atom interface

## Scope and evidence

This is a planning comparison, not a GPAW benchmark or resource release. The exact source geometries were read from the archived extxyz files. No SCF, LAMMPS, MD, or fitting was run. The previously archived resource review documents a 32-GiB cgroup with only about 5.62 GiB headroom at one snapshot, a historical peak within 73,728 bytes of its limit, 15 `memory.max` events, and a 394-atom GPAW initialization that exceeded 32 GiB. There is no matched successful peak-memory measurement for the 217-atom model, and none for M77.

## Geometry and reciprocal-space proxy

| Model | Atoms | Cell volume (Å³) | Proposed mesh | k-point count | Geometry status |
|---|---:|---:|---:|---:|---|
| M77, each of three gaps | 77 | 2041.024552 | 4×4×1 | 16 | Fixed, unrelaxed, 0° shared-strain candidate |
| Archived interface parent | 217 | 5888.439935 | 2×2×1 | 4 | Fixed input; this mesh is also not demonstrated converged |
| M77 optional central-gap convergence check | 77 | 2041.024552 | 5×5×1 | 25 | Proposal only; not run or released |

The cell-volume ratio is `2041.024552 / 5888.439935 = 0.3466155`; the atom-count ratio is `77 / 217 = 0.3548387`. At fixed 500 eV cutoff, a crude plane-wave coefficient/work proxy proportional to `atom_count × cell_volume × number_of_kpoints` gives:

- M77 at 4×4×1 versus 217 atoms at 2×2×1: about **0.492×** the proxy.
- M77 at 5×5×1 versus 217 atoms at 2×2×1: about **0.769×** the proxy.
- At the same k mesh, the corresponding proxy ratio would be about **0.123×**.

These ratios are only scaling indicators. They omit bands and eigensolver workspace, FFT grids, PAW augmentation, parallel distribution, allocator overhead, I/O, checkpointing, and implementation details. They do not predict absolute RAM, elapsed time, or whether a job will fit a particular cgroup. The M77 4×4×1 mesh has more in-plane k points than the archived 217-atom 2×2×1 proposal; it has not been validated as converged.

## Proposed task count and outputs

The package defines three fixed-geometry single-point inputs, one for each prescribed gap (2.16, 2.40, and 2.64 Å), all with the same 4×4×1 mesh and method. A later convergence comparison at the central gap would be an additional, separately authorized SCF at 5×5×1. This task ran zero DFT points. No exact wall-time, disk, or peak-RAM forecast is supported by present measurements.

## Resource decision

```text
M77_RELATIVE_PROXY=LOWER_THAN_217_ATOM_BASELINE
M77_ABSOLUTE_RAM=UNMEASURED
M77_WALLTIME=UNMEASURED
PERSISTENT_CHECKPOINT_CAPACITY=UNVERIFIED
DFT_EXECUTION_RESOURCE_GATE=HOLD
```

Before a later run, require a fresh per-host memory/cgroup and disk preflight, monitored peak RSS/cgroup use, adequate persistent scratch/checkpoint/restore space, four MPI ranks, one OpenMP thread, and a stop threshold that preserves enough reserve for cleanup and archival. If the resource gate fails, stop without reducing cutoff, convergence criteria, or required k-point sampling.
