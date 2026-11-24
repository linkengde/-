# Ag–W TTM reference prototype

This is an Ag–W-specific adaptation of the repository's existing Python finite-volume operator at `models/ag_ti3sic2_qsgs_1to1_finer/ag_ttm_reference.py`. It reuses that solver's CIC atom-to-grid map, zero-flux/open-z finite-volume diffusion, adaptive electron substeps, electron internal-energy tables, and conservative lattice kinetic-energy rescaling. The old Ag–Ti₃SiC₂ model files were not modified. No external solver download was needed.

## Current model

- Ag and W have separate electron energy/temperature fields on a shared mesh.
- Ag/W `Ce` and `G` read the selected Lin 2008 digitized baseline tables. Ag `ke(Te,Tl)` uses the room-anchored Lin Drude closure; W `ke(Te,Tl)` uses the Daraszewicz closure.
- The Gaussian source is mapped only onto Ag atoms in the Newton region (40–200 Å), with `I0=5 eV/(ps Å²)`, `R=30 Å`, `l_skin=20 Å`, center `(100,100) Å`, and surface `z=200 Å`. `S(t)` scales the mapped source and is supplied by the caller.
- Electron diffusion is phase-local. An exploratory `h_e=1 GW m⁻² K⁻¹` face conductance transfers energy equally and oppositely between unlike-phase cell faces. This value uses Au/Pt and Au/Al literature analogues; it is not an Ag/W measurement. The current mesh rule omits sub-cell interface area.
- Electron-to-lattice energy is returned separately for Ag and W and is applied through the reused conservative velocity-rescaling operator.
- The grid spans the Newton layer, with periodic x/y and zero-normal-flux z faces. Phase fractions are recomputed from atom volumes at setup and must sum to no more than one in each cell.

The Ag room-temperature conductivity anchor (429 W m⁻¹ K⁻¹) uses the user-provided bulk total as an electronically dominated proxy. The XTANT 16.56 W m⁻¹ K⁻¹ value remains archived but is not used as the room-temperature amplitude. Ag high-Te conductivity and the interface model remain exploratory.

## Numerical smoke check

Run from the repository root:

```bash
python3 research/ag_w_reference_pack/models/ag_w_qsgs_200A/ttm/run_agw_ttm_grid_smoke.py
```

The smoke check uses a small synthetic two-phase mesh and a deliberately reduced test source. It verifies parameter loading, the 300 K conductivity anchors, Ag-only source mapping, equal/opposite Ag/W interface exchange, and the combined electron-plus-lattice energy ledger, including empty-phase cells in the lattice-temperature arrays. It is not a physical pulse, LAMMPS run, or melt/quench result.

## LAMMPS coupling entry point

`agw_ttm_md_pilot.py` adapts the earlier pure-Ag LAMMPS pilot pattern: EFS `eam/fs`, x/y-periodic and free-z boundaries, fixed/bath/Newton layers, one Ag/W electron grid per phase, and the reused deterministic electron-to-lattice kinetic-energy handoff. By default it preserves the supplied restart without an extra NVT segment, so source/no-source runs can start from the same exact state; an optional shared 300 K equilibration can be requested. It requires explicit pulse FWHM, center, cooling duration, a matching `data_sha256`, and by default a static audit with `md_start_gate: "PASS"`. An explicitly separate `--exploratory` path accepts only a manifest marked `EXPLORATORY_ONLY` that also discloses `STATIC_EFS_FAIL`; it reports a software-link status and must not be interpreted as passing the static or physical model gate.

Example invocation for a future statically accepted structure:

```bash
python research/ag_w_reference_pack/models/ag_w_qsgs_200A/ttm/agw_ttm_md_pilot.py \\
  --data PATH_TO_ACCEPTED.data --audit PATH_TO_STATIC_AUDIT.json \\
  --output PATH_TO_PILOT.json --pulse-fwhm-ps PULSE_FWHM \\
  --pulse-center-ps PULSE_CENTER --cooling-ps COOLING_TIME
```

The run keeps the phase-fraction and atom-to-grid maps fixed after mapping. It is intended only for a short exploratory bridge; substantial melt migration, Ag/W redistribution, and a moving free surface are not modeled. Atoms outside the 40–200 Å material range are excluded from the TTM map, even if they remain in the Newton integration group. Joint mapping moves Ag/W phase volume, Ag source power, and atom-to-grid weights together into nearby cells with shared free capacity; the transfer distance is bounded and reported, and unresolved overfill is rejected. The source is a peak-normalized temporal Gaussian, not normalized to a measured absorbed fluence. The 20–40 Å Langevin bath exchanges energy with the atomistic region; the report records its tallied energy separately from the laser input. Optional `--restart-dump` and `--restart-step` arguments restore an atomistic state after group membership is defined from the template data. The 40–300 Å Newton integration group keeps atoms in the monitor/vacuum region dynamic; only the 40–200 Å material subset couples to the TTM grid.

The Ag/W electronic interface value defaults to a provisional `1 GW m⁻² K⁻¹` from non-Ag/W analogues, and omits sub-cell interface area. The temperature rescaling bridge conserves kinetic energy but does not preserve local momentum or add stochastic phonon fluctuations.

## v6 link preflight

The synthetic grid smoke check passed with a combined energy residual near `1.48e-14 eV`. The v6 finite-size LAMMPS/TTM pilot was attempted from the corrected NVE restart but did not reach a TTM step: the 10 Å atom-to-grid phase fill peaked at 1.31585, and the mapper correctly rejected the overfilled mesh instead of clipping fractions. Coarser tested meshes also remained overfull; only a single-cell whole-domain mapping passed, which cannot resolve phase interfaces or the Gaussian spot and was rejected as unusable. No laser pulse was applied and no matched source/no-source coupled trajectories were completed. See `AGW_TTM_LINK_PREFLIGHT_V6.md` for the mapping scan, environment/debug notes, and retained logs.

## v7 real-structure mapping diagnosis

The v7 diagnosis froze `nve_fixedids_after_nvt_final.dump.gz` and rebuilt Newton membership from the template data-file coordinates before applying the restart. It found 4,490 static Newton-group atoms outside the 40–200 Å TTM material range at the frozen endpoint; 1,875 were already outside in the template and 2,619 moved outside during the NVT/NVE continuation. They are now excluded before CIC, preventing free-surface clipping from concentrating them in the boundary cell. A synthetic boundary check confirms the exclusion.

The new joint Ag/W capacity mapper conserves Ag volume, W volume, Ag source power, and atom-to-grid weights while rejecting unresolved overfill. Its synthetic transfer check and the separate synthetic TTM energy-ledger check pass. However, the actual 5, 10, and 20 Å maps all fail the 20 Å local capacity limit. A maximum-flow feasibility diagnostic also finds 6,966, 15,711, and 19,186 Å³ of overflow that cannot reach receiver capacity within that bound, respectively. Therefore no real-structure TTM step and no source/no-source LAMMPS/TTM pair has been run. See [`AGW_TTM_MAPPING_V7.md`](AGW_TTM_MAPPING_V7.md) and `agw_ttm_mapping_v7_report.json` for the counts and full scan.

The source quadrature on the frozen coordinates is 9,907.477809 eV/ps at `S(t)=1`; the integral is consistent across the tested mesh spacings, but no heat pulse was applied. The v5 Dai EFS static audit remains failed, and the interface electronic conductance is still an exploratory analogue value.

For an end-to-end LAMMPS mapping gate without starting a coupled trajectory, `agw_ttm_md_pilot.py --mapping-only` performs the setup and mapping, writes a PASS/FAIL JSON, and exits before constructing the electron solver or integrating a TTM step. On the current frozen restart at 10 Å, this gate returned FAIL with 107.925388 Å³ unplaced at cell 176, recorded zero TTM steps, and applied no laser source; see `agw_ttm_mapping_v7_lammps_preflight.json` and its adjacent LAMMPS log. This entry point can be rerun after a future geometry or mapping change.
