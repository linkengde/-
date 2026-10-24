# B1 mesh execution pack — 2026-10-24

This folder turns the phase-specific mesh plan into twelve exact endpoint records for six same-geometry comparisons. It is a **preflight-only pack**: it has no GPAW calculation runner, every launch and integration flag is false, and verification never writes labels.

The GPAW/ASE stack is already installed, so no software download was needed. Verified versions are GPAW 26.7.0, ASE 3.29.0, and PAW data 1.2.1. `gpaw info` reports MPI and ScaLAPACK enabled. The Ag, Ti, Si, and C PAW setup SHA256 values match the common method proposal. This validates availability, not a live four-rank calculation.

The manifest covers Ag32, Ag64, Bai Ti₃SiC₂ 72/96, and liquid Ag 1250/1600 K endpoints. Every label points to one frozen geometry hash; each pair differs only in k-point mesh. Energy target is `free_energy_eV`, with the native extrapolated energy retained separately. Liquid 3×3×3 jobs remain conditional and are added only if a Γ/2×2×2 pair misses a threshold.

Run the safe, read-only preflight from the repository root:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/coordination/reports/window-c/b1_kmesh_execution_pack_20261024/verify_pack.py
```

This script checks all input hashes and frozen candidate roles, verifies the twelve rows form six matched pairs, rejects any existing DFT label, and verifies the local GPAW/ASE/PAW installation. It does not start MPI or GPAW. Expected status is `PASS_PACK_VALIDATED_LAUNCH_DISABLED`.

The `σ=0.10 eV` method remains provisional until Window A publishes and verifies the registered V20 `σ=0.05` and `σ=0.20` controls. The last task-card snapshot is stale and there is no live owner telemetry in this checkout. Once the owner gate is complete, freeze the final method, register these controls with an owner and resource limit, then add the execution runner. Until then the pack remains disabled.
