# Ag–W V9 hybrid reconstruction audit

**Overall status: FAIL — no V9 TTM-ready candidate was produced.** V9 preserved the full atom structure and completed the requested mesh, cold-baseline, and LAMMPS mapping checks. The results are saved as an auditable failed milestone; no atom deletion, box expansion, atom-volume change, or capacity-limit relaxation was used.

## Scope and frozen inputs

V9 combines V7's complete V7-derived atom structure with V8's grid-phase search and joint Ag/W mapping goal. The following inputs were frozen:

- Template: `interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz`, SHA-256 `29bf0eccea5516f1e51e7454876d7f66a99fd3b4fcdb3905b2e0e27993d3bde8`.
- Phase-scan restart: V6/V7 cold candidate at timestep 4200, SHA-256 `1814c0f9c26831b66ff85474ed1612ffd1ed2a76913e7493d96b1ca9f6b5abcc`.
- Full structure: 482,601 atoms (Ag 304,046; W 178,555), Ag/W mass ratio 0.999127.
- Box and boundary: 200×200×300 Å, p p f.
- Fixed IDs: the original 46,117 atoms in the bottom 0–20 Å layer.
- TTM atomic-volume constants: unchanged from the existing Ag–W reference solver.

The phase scan used the V6/V7 frozen restart first. A second 10 Å phase scan and LAMMPS preflight used the later V9 cold-continuation restart. The full structure, box, IDs, and total species counts were preserved in both.

## Stage results

| Stage | Result | Evidence |
|---|---|---|
| Preserve full structure and original box | **PASS** | Input SHA, atom IDs/types, atom count, composition, and box checked |
| Scan mesh spacing and x/y phase | **COMPLETE; no feasible mapping** | 84 phase combinations scanned for 5/10/20 Å; 10/20 Å phases use the shared direct-flow routine. For 5 Å, all four raw CIC phases were scanned and the lowest-overflow phase was sent to the larger flow graph. |
| Share capacity diagnostic and production algorithm | **PASS as implementation** | `joint_capacity_maxflow` is called by both `build_agw_grid_fields_joint` and the V7 diagnostic wrapper. The LAMMPS mapping-only path confirms the production result. |
| 300 K cold NVT/NVE continuation | **FAIL as a balanced baseline** | Two 0.5 ps NVT plus 0.05 ps NVE segments completed. Global temperature and NVE energy behavior are stable, but Ag/W temperatures separate and potential energy continues to fall. |
| Actual 10 Å LAMMPS mapping-only preflight | **FAIL** | Best post-cold phase (5,5) Å leaves 6,533.52831 Å³ unassigned within 20 Å. Zero TTM steps; no laser source. |
| Source/no-source TTM pair | **NOT RUN** | Blocked by the failed mapping and cold-baseline gates. |

## Full-structure phase scan

The x/y phase increment was 2.5 Å. Grid spacings were 5, 10, and 20 Å; transfer distance was fixed at 20 Å. Each donor can transfer only to a cell that had initial free capacity. No transferred material can become a second donor.

| Grid | Direct-flow candidate phase | Raw overflow demand (Å³) | Unassigned (Å³) | Raw peak combined fill | Laser-profile note |
|---:|---:|---:|---:|---:|---|
| 5 Å | (0, 2.5) Å | 99,952.695 | 6,994.000 | 1.309 | Best raw-overflow phase; still fails. The other 5 Å phases were retained in the scan, but not sent through the costly flow graph. |
| 10 Å | (0, 2.5) Å | 48,476.070 | 15,655.717 | 1.119 | Best direct-flow phase among all 16. Source radial RMS 29.40 Å; depth mean 21.52 Å. |
| 20 Å | (17.5, 17.5) Å | 27,921.563 | 18,774.187 | 1.062 | Best direct-flow phase among all 64, but too coarse: R=30 Å spans 1.5 cells and the first 0–10 Å radial bin receives zero source power. |

Global grid volume is sufficient, but local free volume within the allowed transfer distance is not. These values describe the chosen atomic-volume/CIC representation and do not alone prove that the atoms are physically compressed.

The V9 capacity solver uses integer maximum-flow capacities. It selects a precision that keeps the SciPy flow values within the backend's integer range; it then places any fractional rounding tail only into original receiver cells with actual remaining capacity. The quantization bound is reported as a diagnostic and is **not** treated as extra cell capacity. A phase passes only when residual donor fill is within the same per-cell numerical tolerance enforced by the final mapped-fraction check.

## Post-cold 10 Å phase scan

The final V9 restart after the second NVE segment was rescanned across all 16 phases using the same shared production capacity solver:

- Best phase: (5,5) Å.
- Raw overflow demand: 30,559.16835 Å³.
- Unassigned within 20 Å: 6,533.52831 Å³.
- Peak raw fill: 1.10326.
- Mapped Newton-layer atoms: Ag 235,272; W 145,383.
- 7,619 atoms in the original static Newton group were outside the 40–200 Å TTM mapping region and excluded instead of being clipped into boundary cells.
- Ag source quadrature: 10,015.18054 eV/ps. This is a mapping diagnostic only; it was not applied.
- The best-phase source RMS radius is 29.22 Å. At 20 Å radial binning, the first 0–10 Å bin gets 4.52% of power; a coarse 20 Å grid instead gave zero in that bin.
- Full model mass ratio remains 0.999127. The mapped 40–200 Å slab has a different local species ratio because of the spatial composition and excluded atoms; it must not be confused with the whole-model 1:1 target.

Actual LAMMPS mapping-only preflights were run at the V8-recommended (0,2.5) Å phase and the best post-cold (5,5) Å phase. They returned the same shared-solver failures, respectively 6,703.6155 Å³ and 6,533.52831 Å³ unassigned. Both ran the LAMMPS setup/force screen only. Each report records `ttm_steps_executed: 0` and `laser_source_applied: false`.

The final LAMMPS `run 0` force screen was max 3.9031 eV/Å and RMS 0.8579 eV/Å. It is diagnostic only; the original static EFS audit remains `STATIC_EFS_FAIL`.

## Cold baseline

The frozen V7-derived restart was continued with a single 300 K NVT thermostat on all mobile atoms. The original bottom 46,117 fixed IDs were kept. Each segment used 0.25 fs steps, 0.5 ps NVT, then 0.05 ps NVE. No phase-specific thermostat, TTM, or laser source was used.

| Segment | NVT mobile / Ag / W end temperature (K) | Last 0.2 ps mobile ΔT (K) | NVT potential-energy change (eV) | NVE energy drift (eV) |
|---|---:|---:|---:|---:|
| First | 302.77 / 295.20 / 315.18 | +0.40 | −576.76 | −0.00036 |
| Follow-up | 303.54 / 293.21 / 320.47 | −1.31 | −902.95 | −0.00001 |

All 482,601 atoms and IDs were retained; fixed-ID displacement was zero. The final mobile maximum displacement over the follow-up segment was 2.575 Å, with 3,440 mobile atoms moving more than 1 Å. In the subsequent 0.05 ps NVE check, Ag/W ended at 293.73/323.02 K. The stable global temperature and excellent NVE energy conservation do not remove the growing 27.26 K Ag/W temperature gap or the continuing potential-energy relaxation. The cold-baseline gate therefore remains **FAIL** for a common 300 K Ag–W starting state.

## Interpretation and next gate

V9 demonstrates that grid-phase optimization and a direct joint capacity solver improve diagnostics over the previous greedy/cascading mismatch, and that NVT/NVE can run without atom loss. It does **not** produce a usable TTM start. The direct capacity scan still cannot place all estimated phase volumes locally, and the cold state remains phase-temperature imbalanced. No source/no-source pair, laser heating, melting, or solidification run was started.

Three-dimensional Ag/W percolation was not reclassified after the 1 ps NVT coordinate evolution; this audit does not claim that the dynamically evolved connectivity gate passed.

The next technical direction is to inspect the local density representation using a physical geometry-based measure (for example, validated Voronoi volumes) and/or reconstruct interface geometry. Any replacement must preserve total Ag/W volume and direct source power while checking actual per-cell capacity. Do not label geometry-mask fractions as atom-derived TTM fractions without a separate consistency check. Re-run the cold baseline and actual LAMMPS mapping gate after any new structure or mapping method.

## Reproduction

Run from the repository root with the Ag–W environment:

```bash
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v9_hybrid/scan_v9_resolution_phase.py
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v9_hybrid/verify_v9_shared_flow.py
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v9_hybrid/analyze_v9_cold_baseline.py
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v9_hybrid/analyze_v9_cold_followup.py
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v9_hybrid/scan_v9_postcold_10A.py
```

The LAMMPS inputs and complete mapping-only JSON/log records are in this directory. They require the existing `agw_ttm_md_pilot.py` entry point and the V9 exploratory manifest. The manifest explicitly retains the failed static EFS audit and does not authorize TTM integration.

The actual cold continuation and final-state mapping gate were run with:

```bash
cd research/ag_w_reference_pack/models/ag_w_qsgs_200A
OMP_NUM_THREADS=4 LD_LIBRARY_PATH=/workspace/.venvs/ag-qsgs/lib:${LD_LIBRARY_PATH:-} \
  /workspace/.venvs/ag-qsgs/bin/lmp -pk omp 4 \
  -in interface_reconstruction_v9_hybrid/in.v9_cold_nvt_nve \
  -log interface_reconstruction_v9_hybrid/v9_cold_nvt_nve.log
OMP_NUM_THREADS=4 LD_LIBRARY_PATH=/workspace/.venvs/ag-qsgs/lib:${LD_LIBRARY_PATH:-} \
  /workspace/.venvs/ag-qsgs/bin/lmp -pk omp 4 \
  -in interface_reconstruction_v9_hybrid/in.v9_cold_followup \
  -log interface_reconstruction_v9_hybrid/v9_cold_followup.log
cd ../../../../
LD_LIBRARY_PATH=/workspace/.venvs/ag-qsgs/lib:${LD_LIBRARY_PATH:-} \
  /workspace/.venvs/ag-qsgs/bin/python research/ag_w_reference_pack/models/ag_w_qsgs_200A/ttm/agw_ttm_md_pilot.py \
  --data research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz \
  --audit research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v9_hybrid/v9_exploratory_state_manifest.json \
  --restart-dump research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v9_hybrid/v9_nve2_end.dump.gz \
  --restart-step 2200 --output /tmp/agw_v9_mapping.json \
  --pulse-fwhm-ps 0.01 --pulse-center-ps 0.02 --cooling-ps 0.02 \
  --mapping-only --exploratory --grid-A 10 \
  --grid-offset-x-A 5 --grid-offset-y-A 5 --max-map-transfer-A 20
```

The last command only performs LAMMPS setup and a mapping preflight. Its expected exit is a recorded mapping FAIL; it performs no TTM time step and applies no source.
