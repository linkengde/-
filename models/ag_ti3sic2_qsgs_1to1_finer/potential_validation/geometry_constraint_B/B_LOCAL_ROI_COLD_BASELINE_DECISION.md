# Local-ROI Model B cold-baseline decision

Date: 2026-10-25. The whole-layer Model B baseline had failed. A single local-ROI variant was then tested to isolate the source neighborhood from distant Ag rearrangement. This was a new boundary-condition definition, not another relaxation cycle on the whole layer.

## Model and test

- Reused the 916,534-atom relaxed-parent QSGS geometry and archived Ag-only FIRE-100 coordinates.
- Kept Ti/Si/C rigid and fixed all Ag outside a cylinder of radius 75 Å centered under the heat source. Only 33,963 Ag atoms inside the cylinder and above z=20 Å were mobile.
- Initialized mobile Ag at 300 K and ran a no-source NVE test for 500 steps at 0.5 fs (0.25 ps).
- The operational source core was defined as Ag within lateral radius 30 Å and surface depth 0–60 Å; it contains 1,275 atoms.
- Reproduction from the finer-model directory:

```bash
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp \
  -in potential_validation/geometry_constraint_B/in.local_roi_cold_nve_0p25ps \
  -log potential_validation/geometry_constraint_B/log.local_roi_cold_nve_0p25ps
python potential_validation/geometry_constraint_B/analyze_local_roi_cold_nve.py
```

## Results

| Metric | Start | 0.25 ps | Change |
|---|---:|---:|---:|
| Mobile-ROI Ag temperature | 300.00 K | 224.97 K | −75.03 K |
| Potential energy | −537,983.48 eV | −537,654.11 eV | +329.37 eV |
| Kinetic energy | 1,316.98 eV | 987.62 eV | −329.36 eV |
| Total energy | −536,666.50 eV | −536,666.49 eV | +0.01 eV |
| Maximum mobile-Ag force | 1.348 eV/Å | 3.542 eV/Å | increased |

The 33,963-atom mobile ROI moved 0.283 Å RMS (p90 0.425 Å; p99 0.927 Å). The 1,275-atom source core moved 0.380 Å RMS (p90 0.590 Å; p99 1.050 Å), while its force p95/max increased from 0.757/1.090 to 0.959/1.897 eV/Å.

The small total-energy drift shows that the integration is behaving well. The large opposing potential/kinetic changes, falling temperature, and rising forces show that the local structure is still adjusting without any applied heat. Freezing distant Ag therefore did not produce a stationary cold baseline for the target region.

## Decision

**No-go for adding TTM or a heat source to this local-ROI B variant.** The current whole-layer and local-ROI geometry-only preparations both fail their cold baseline checks. This bounded test completes the current B attempt; do not extend it with further FIRE/NVT cycles or interpret a heated trajectory from this initial state as melting.

The target core has smaller rearrangement than the full ROI, but it still does not start from a demonstrated stable 300 K baseline. Any next attempt needs a materially different initial-state construction or a better-supported Ag–scaffold model, with its own finite cold-state gate. The rigid scaffold/WCA model remains useful only as an archived exploratory geometry test; it has not formed or validated a melt pool.
