# Window A B1 phase-specific mesh task — 2026-10-26

This is a separately versioned Window A task for the six fixed-geometry k-point convergence pairs proposed in the B1 phase-specific mesh plan. The Bai candidate package checksum inventory was corrected for three stale entries; the complete inventory now verifies. The original twelve-endpoint preflight also passes.

## Frozen scope and method

The task runs the twelve endpoints in `task_manifest.json` sequentially, four MPI ranks at a time. Each pair uses identical input bytes, cell, atom order, positions, PBC, XC, plane-wave cutoff, Fermi smearing, mixer, and SCF thresholds; only the k-point mesh changes. The common task method is GPAW 26.7.0 / ASE 3.29.0 / PBE / PAW 1.2.1 / PW 500 eV / sigma 0.10 eV, with the force-consistent `free_energy` target retained as `free_energy_eV` and GPAW's native extrapolated energy recorded separately.

This freezes the method for these mesh comparisons only. It does not claim a universal correction between bulk and slab settings, change any approved dataset, or finalize the method for a successor model.

The six mandatory pairs are:

| Geometry | Coarse mesh | Fine mesh |
|---|---:|---:|
| Ag32, seed 9073 | 6×6×6 | 8×8×8 |
| Ag64, seed 9073 | 5×5×5 | 7×7×7 |
| Bai Ti₃SiC₂ 72 atoms, seed 9073 | 5×8×2 | 7×10×2 |
| Bai Ti₃SiC₂ 96 atoms, seed 9073 | 8×8×2 | 10×10×2 |
| Liquid Ag, 1250 K, seed 501 | Γ | 2×2×2 |
| Liquid Ag, 1600 K, seed 902 | Γ | 2×2×2 |

The two 3×3×3 liquid endpoints are conditionally registered and may run only if the corresponding Γ→2×2×2 comparison exceeds a threshold. The total cap is fourteen endpoints.

## Resources and acceptance

- One endpoint at a time; four MPI ranks per endpoint; no GPU; no automatic retries.
- Compare each pair before proceeding. Record the free-energy difference per atom, force-vector RMS difference, and maximum atom-force-vector difference.
- Thresholds are 2 meV/atom, 0.01 eV/Å RMS, and 0.03 eV/Å maximum, respectively.
- Archive compact result extxyz, summary, progress, GPAW log, verification receipt, and SHA256 inventory. Keep large GPAW state files local and remove them only after the compact archive passes verification.
- Outputs remain numerical method controls with `TRAINING_CANDIDATE_ONLY` input lineage. They are not integrated into V19 or any other dataset and cannot be used for training or validation without a later, separate review.

## Ag64 7×7×7 OOM recovery

The first Ag64 7×7×7 attempt was killed with exit code 137 during GPAW's LCAO-to-PW initialization, before SCF began. The cgroup recorded an OOM kill at its 32 GiB memory limit. Its logs and progress record are preserved under `failed_attempts/`.

A manually registered recovery used GPAW's random plane-wave initial wavefunctions for this endpoint only, while retaining four MPI ranks and the physical method settings. It was also OOM-killed before SCF at the 32 GiB limit. Both failures are preserved in `failed_attempts/`; the queue is stopped pending a reviewed memory/parallelization strategy. Automatic retries remain disabled, and the task manifest launch flag is disabled until resource review is complete.

## Execution

Run the task verifier and existing source-pack verifier from the repository root before launch:

```bash
python research/mace-v12-transfer/coordination/reports/window-a/b1_phase_specific_mesh_task_20261026/verify_task.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/coordination/reports/window-c/b1_kmesh_execution_pack_20261024/verify_pack.py
```

The registered queue is `run_queue.py`. It takes the shared four-rank DFT lock, launches labels in the listed pair order, verifies each completed archive, compares each pair, and conditionally runs a liquid 3×3×3 endpoint when required. The queue never edits dataset files or starts model training.
