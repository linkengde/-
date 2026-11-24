# Ag–W TTM mapping v7 diagnosis

## Decision

**Real-structure mapping is still rejected at the configured 20 Å maximum cell-transfer distance. Do not run the source/no-source LAMMPS/TTM pair yet. No TTM step or laser source was applied in this diagnosis.**

This is a mapping-capacity result for the frozen v6 restart, not evidence that the atomic material is physically compressed by the reported cell-fill peaks. Global phase volume fits within the total grid capacity, but enough free space is not reachable from every overfilled cell within the 20 Å limit.

## Frozen inputs and region membership

- Template data: `interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz` (`SHA-256: 29bf0eccea5516f1e51e7454876d7f66a99fd3b4fcdb3905b2e0e27993d3bde8`).
- Frozen restart: `interface_reconstruction_v6_exploratory/nve_fixedids_after_nvt_final.dump.gz`, timestep 4200 (`SHA-256: 1814c0f9c26831b66ff85474ed1612ffd1ed2a76913e7493d96b1ca9f6b5abcc`).
- LAMMPS group membership is established from template coordinates before `read_dump` replaces coordinates. The diagnostic now reproduces that static Newton group from the same data file; it does not infer group membership from evolved z coordinates.
- Static Newton group: 388,274 atoms (240,826 Ag and 147,448 W). The template itself already placed 1,875 of them outside the TTM mapping interval. At the frozen restart, 383,784 remain in the interval 40–200 Å (237,404 Ag and 146,380 W); 4,490 are outside it (3,422 Ag and 1,068 W), including 2,619 atoms that were inside the interval in the template and moved outside by the NVT/NVE continuation.

The fixed-ID Newton group had 4,474 atoms outside 40–200 Å at the NVT endpoint (2,603 newly outside the template interval) and 4,490 at the subsequent NVE endpoint (2,619 newly outside). At the frozen restart, 4,384 atoms are above z=200 Å. The former z-clipped CIC code would have placed all such above-surface atoms in the top mesh layer; those atoms are now filtered before CIC. If the legacy source were evaluated without that filter, the above-surface Ag atoms alone would contribute about 1,040.41 eV/ps at `S(t)=1` under the fixed z=200 Å source plane. That value is diagnostic only; it was not injected.

## Raw CIC and capacity scan

The raw CIC fill includes atomic-volume proxies of 4.09³/4 Å³ for Ag and 3.16³/2 Å³ for W. The joint mapper uses periodic x/y neighbors and open z neighbors. The 20 Å limit is an exploratory spatial-locality screen for this audit, not a measured material length or calibrated mapping parameter. The maximum-flow diagnostic is an independent feasibility bound with 0.001 Å³ integer capacity resolution; it checks whether donor overflow can be assigned to real receiver free volume within that limit, without relying on the production mapper's greedy donor order.

| Grid spacing | Raw maximum fill (Ag / W / combined) | Raw combined overfilled cells | Overflow demand | Global free volume | Max-flow assignable | Unassignable within 20 Å |
|---:|---:|---:|---:|---:|---:|---:|
| 5 Å | 1.2999 / 1.1613 / 1.2999 | 24,997 | 100,897.07 Å³ | 130,751.01 Å³ | 93,930.61 Å³ | 6,966.46 Å³ |
| 10 Å | 1.1202 / 1.0685 / 1.1202 | 2,899 | 48,255.04 Å³ | 78,108.97 Å³ | 32,543.59 Å³ | 15,711.45 Å³ |
| 20 Å | 1.0625 / 1.0113 / 1.0636 | 276 | 28,114.16 Å³ | 57,968.10 Å³ | 8,928.17 Å³ | 19,185.99 Å³ |

All three grids fail the same local-capacity condition even under the maximum-flow feasibility check. This makes the 10 Å failure unrelated to donor ordering in the greedy implementation. Independent per-phase redistribution is also invalid: its combined peak fill remains 1.7428, 1.3160, and 1.0858 for 5, 10, and 20 Å grids respectively.

The Ag source quadrature on the mapped restart coordinates is 9,907.477809 eV/ps at `S(t)=1`; the atom-wise integral and raw CIC integral agree to floating-point precision on all three meshes. This verifies the source integral only. It is not a delivered heat pulse.

## Algorithm and checks

`agw_ttm_reference.py` now has a joint Ag/W capacity mapper. It transfers the same fractional package of Ag volume, W volume, Ag source power, and atom-to-grid weights from each donor cell into nearby cells with available capacity. It records transfer distances, preserves phase volumes and source power, and rejects unresolved overfill rather than clipping or renormalizing fractions.

The synthetic joint-transfer check passed: combined fill stayed at or below one; Ag and W volumes and Ag source power were conserved; atom-map rows remained normalized; and observed transfer stayed within the configured distance. The synthetic free-surface check passed: an atom above z=200 Å was excluded before CIC and contributed no source. The synthetic-grid TTM energy-ledger smoke check also passed with a residual of `1.48e-14 eV`; it does not validate the real structure map.

The LAMMPS runner was corrected to capture the Newton-group mask from the template data before loading restart coordinates, matching the static group membership used by LAMMPS. Mapping then uses only those group members whose current coordinates remain in 40–200 Å. Atoms outside that interval stay in the Newton integration group but are excluded from the material TTM mesh.

An end-to-end LAMMPS `--mapping-only` preflight was run for the 10 Å mesh on the same data/restart hashes. LAMMPS setup and `run 0` completed; mapping stopped with 107.925388 Å³ of overflow remaining at cell 176 under the 20 Å limit. The output records `AGW_MAPPING_PREFLIGHT_FAIL`, `ttm_steps_executed: 0`, and `laser_source_applied: false` in `agw_ttm_mapping_v7_lammps_preflight.json`.

## Reproduction and next gate

From the repository root, using the configured `ag-qsgs` Python environment:

```bash
/workspace/.venvs/ag-qsgs/bin/python research/ag_w_reference_pack/models/ag_w_qsgs_200A/ttm/diagnose_agw_ttm_mapping_v7.py
/workspace/.venvs/ag-qsgs/bin/python research/ag_w_reference_pack/models/ag_w_qsgs_200A/ttm/run_agw_ttm_grid_smoke.py
```

The LAMMPS integrated mapping gate used for this restart can be repeated from the repository root with:

```bash
LD_LIBRARY_PATH=/workspace/.venvs/ag-qsgs/lib /workspace/.venvs/ag-qsgs/bin/python research/ag_w_reference_pack/models/ag_w_qsgs_200A/ttm/agw_ttm_md_pilot.py \
  --data research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz \
  --audit research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v6_exploratory/exploratory_ttm_state_manifest.json \
  --exploratory --restart-dump research/ag_w_reference_pack/models/ag_w_qsgs_200A/interface_reconstruction_v6_exploratory/nve_fixedids_after_nvt_final.dump.gz \
  --restart-step 4200 --output research/ag_w_reference_pack/models/ag_w_qsgs_200A/ttm/agw_ttm_mapping_v7_lammps_preflight.json \
  --pulse-fwhm-ps 0.01 --pulse-center-ps 0.02 --cooling-ps 0.02 --grid-A 10 --max-map-transfer-A 20 --mapping-only
```

`--mapping-only` requires zero equilibration and exits before constructing the electron solver or integrating a coupled step.

Machine-readable results are in `agw_ttm_mapping_v7_report.json`. The next step is to revise the atom-to-grid representation or rebuild the density field while preserving the capacity guard, then repeat the 5/10/20 Å conservation and resolution checks. Increasing the transfer radius without a spatial-smearing assessment is not treated as a pass. The Dai EFS static audit remains `STATIC_EFS_FAIL`; the interface conductance remains an exploratory analogue value; and the phase fractions/map remain fixed during a short prototype run.
