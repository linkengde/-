# Ag-W QSGS candidate, 200 Å material region

This is a **new, unrelaxed Ag-W geometry candidate**. It does not modify the Ag-Ti₃SiC₂ model. The QSGS-grown W phase is the backbone; Ag fills the complementary voxel phase.

## Dimensions and layer groups

| Item | Value |
|---|---:|
| Material region | 200 × 200 × 200 Å |
| Simulation box | 200 × 200 × 300 Å |
| Boundary | x/y periodic, z nonperiodic (`p p f`) |
| Vacuum | z = 200–300 Å |
| Fixed layer | z = 0–20 Å |
| Thermostat layer | z = 20–40 Å |
| Newton layer | z = 40–200 Å |
| Sputter monitor | z = 200–250 Å |
| Evaporation monitor | z = 250–300 Å |

The two monitor regions are empty at the initial state. They are vacuum-space bins for later trajectories.

## Four-parameter QSGS and composition

The builder defines its four QSGS controls explicitly: core seed probability 0.0003, base growth probability 0.145, maximum 10 generations, and target W voxel fraction 0.3545. W grows through 26 neighboring directions with probability `p_grow / r`; x/y are periodic and z is open. Ag occupies the complementary mask.

The resulting 80 × 80 × 80 voxel morphology has 2.5 Å voxels and 35.45 vol.% W / 64.55 vol.% Ag before atomic lattice mapping. This is close to the density estimate of 35.27 vol.% W / 64.73 vol.% Ag for equal masses.

| Item | W | Ag |
|---|---:|---:|
| Lattice | BCC, a = 3.16 Å | FCC, a = 4.09 Å |
| Atom count | 178,435 | 304,036 |
| Atom fraction | 36.9836% | 63.0164% |
| In-plane strain | +0.4621% | −0.2046% |

The Ag/W mass ratio is 0.999766 (relative error 0.0234%). Connectivity was checked on the atomized structures using 2.95 Å for W and 3.15 Å for Ag: both largest components have xy winding rank 2 and span the material thickness. The largest component contains 98.59% of W atoms and 99.82% of Ag atoms; small disconnected fragments remain and are counted in the report.

## Files and checks

- `ag_w_qsgs.data.gz`: LAMMPS atomic data, type 1 = W and type 2 = Ag.
- `ag_w_qsgs.extxyz.gz`: coordinates with `phase`, `layer`, `fixed_layer`, and `thermostat_layer` arrays.
- `qsgs_morphology.npz`: W/Ag voxel masks and the QSGS controls.
- `geometry_report.json`: counts, mass ratio, strains, connectivity, layers, and caveats.
- `in.geometry_read_check`: LAMMPS verifies the data file and layer counts without integrating atoms.

Regenerate into a fresh output directory with:

```bash
python research/ag_w_reference_pack/scripts/build_ag_w_qsgs.py --output /path/to/new/ag_w_qsgs_200A
```

**Do not start MD from the raw data file yet.** Independent FCC/BCC lattices create an initial minimum Ag–W distance of 0.08 Å. An EFS static-minimization preflight reduced the maximum force but stopped with 2.17 eV/Å residual force and a 0.37 Å closest Ag–W pair; it did not produce an acceptable starting state. See `overlap_relaxation/preflight_report.md`. No 300 K trajectory or TTM heat pulse has been run.

A follow-up repair audit tried soft-core clearing, same-phase neighbor springs, and a constrained Dai EFS relaxation. Each route still left severe short contacts or high residual forces, so **no repaired geometry passed static validation**. The original QSGS file remains unchanged. See [`overlap_repair_v1/OVERLAP_REPAIR_AUDIT.md`](overlap_repair_v1/OVERLAP_REPAIR_AUDIT.md) for measurements and rejected methods.

### Interface reconstruction v2 status

The next reconstruction audit tested mass-balanced conflict-graph site pruning followed by unconstrained Dai EFS minimization, full-atom pair projection, sequential pair projection, and same-phase lattice springs followed by EFS. **All routes were rejected.** The pruned candidate removed 49,669 atoms (10.3%) and the 2,700-step EFS relaxation remained unconverged with 1.969 eV/Å maximum force norm. Projection transferred short contacts into same-phase lattices; spring-assisted EFS retained a 0.37 Å Ag–W contact and 11.77 eV/Å maximum force component. The original QSGS input is unchanged.

No full-size QSGS candidate has passed static acceptance. **Do not start 300 K MD, TTM, or laser heating from the raw structure.** The old 9,784.984 eV/ps Ag-only source calculation was made on the raw, overlapping coordinates and is not valid for a future reconstructed structure; recalculate it only after an accepted geometry exists. See [`interface_reconstruction_v2/RECONSTRUCTION_V2_AUDIT.md`](interface_reconstruction_v2/RECONSTRUCTION_V2_AUDIT.md) and its scripts, logs, and machine-readable reports for full evidence.

### Interface reconstruction v3 local slab screen

V3 searched low-strain coincidence cells, then relaxed one-interface seven-layer slabs. It completed 360 runs over all distinct layer-sequence offsets, the preregistered translations, and ten initial gaps. Every scan converged below the 0.03 eV/Å force limit, retained its atoms, and passed the same-phase collapse check. The DFT energy comparison did not pass, most clearly for Ag(111)/W(110). The project accepts the Dai potential for **qualitative graduate-level exploration**, with its bias stated; this does not make the DFT screen a pass or validate quantitative interface energies. See [`interface_reconstruction_v3/LOCAL_INTERFACE_V3_AUDIT.md`](interface_reconstruction_v3/LOCAL_INTERFACE_V3_AUDIT.md), `interface_reconstruction_v3/local_interface_v3_report.json`, and `interface_reconstruction_v3/termination_runs/`.

This local result does not repair the full QSGS atom arrangement. The original 482,471-atom model still has a 0.08 Å cross-phase contact, and no v1/v2 full-size repair passed. A new 3D candidate must preserve atom density, composition, and bicontinuity and pass static checks before any 300 K MD. Ag/W electronic-property tables are now available as candidate sources, but a TTM implementation still needs explicit species/interface closures and pulse conditions.

### V3 exploratory full-size geometry attempt

An all-atom geometric projection kept all 482,471 atoms and the 1:1 mass ratio, but an unconstrained Dai EFS minimization of the resulting candidate reached its 10,000-step limit with 0.650 eV/Å maximum mobile force. It produced a 2.099 Å Ag–Ag minimum, a 2.267 Å Ag–W minimum, and 1,604 atoms above the z=200 Å material surface. The RMS force was low, but the collapsed contacts and surface displacement fail the practical static-start gate. **Do not use this candidate for 300 K MD or TTM yet.** The audit, candidate, input, logs, and analysis scripts are retained in [`interface_reconstruction_v3/qsgs_candidate/FULL_SIZE_STATIC_AUDIT.md`](interface_reconstruction_v3/qsgs_candidate/FULL_SIZE_STATIC_AUDIT.md).

TTM remains the intended heating method after the structural and thermal-model gates are passed. The stock `fix ttm/mod` does not support this free-surface boundary. Ag/W XTANT-3 `Ce`, `G`, and `K` source tables have now been archived, but the W coupling models differ substantially, the supplied conductivity tables depend on `Te` only, and no Ag/W electronic interface term is available. Pulse duration/envelope and absorbed fluence are still unspecified. The full-size geometry has not passed static acceptance, so no MD or TTM trajectory has been run. See `ttm_source_preflight/TTM_SOURCE_PREFLIGHT.md`, the project-level `../../TTM_HEATING_READINESS.md`, and `../../ttm_data/XTANT3_AG_W_DATA_AUDIT_20261124.md` for scope and limits.

### Interface reconstruction v5

V5 regularized the original QSGS mask at fixed W voxel volume and regenerated separate W-BCC and Ag-FCC atom populations. The best geometry-only candidate met the composition, estimated-density, and bicontinuity checks, but seven W-W distances remained below the frozen screen. Spring-free Dai EFS relaxation with the fixed bottom layer stopped at its 10,000-step cap. The final structure had 10,302 W-W, 681 Ag-Ag, and 797 Ag-W pairs below the frozen screens, 1,871 atoms above z=200 Å, and 0.381 eV/Å maximum moving-atom force. **V5 fails strict static acceptance and the exploratory start gate; no 300 K MD or TTM was run.** The TTM solver is unchanged. See [`interface_reconstruction_v5/RECONSTRUCTION_V5_AUDIT.md`](interface_reconstruction_v5/RECONSTRUCTION_V5_AUDIT.md), `precondition_trial_summary.json`, and the retained LAMMPS inputs and logs.

### Exploratory short NVE pilot

Following the revised tiered criteria, a separate v6 pilot used the v5 finite Dai-EFS endpoint, fixed only z=0–20 Å, and ran 25 fs of 300 K no-source NVE at 0.25 fs. It completed with unchanged atom count and energy drift below the log's 0.1 eV display resolution, so the level 1 program screen passed. The cold baseline did not: the mobile temperature fell to 197.5 K as potential energy rose by 5.78 keV, indicating ongoing structural relaxation. **No TTM source, laser, or melt/cool interpretation was applied.** The v5 regularization modifies 8.41% of the QSGS voxels while keeping the source files intact; this candidate is not an exact voxel copy. See [`interface_reconstruction_v6_exploratory/README.md`](interface_reconstruction_v6_exploratory/README.md) and `nve_pilot_metrics.json`.

A follow-on 0.5 ps NVT at 300 K and 0.05 ps NVE check also completed. The system stayed near 197 K through NVT instead of reaching the 300 K target; the thermostat added about 4.44 keV while most net energy accumulated as potential energy. The short NVE segment conserved energy to 0.003 eV near 195 K. The cold-baseline stage therefore remains incomplete, and no TTM or laser run was started. See `interface_reconstruction_v6_exploratory/nvt_followup_metrics.json` and its reproduction input/log.

The corrected fixed-ID continuation then ran for 1.0 ps NVT plus 0.05 ps NVE from the preceding endpoint. It restored the original 46,117 fixed IDs and retained all atoms. Mobile temperature reached 293.79 K; the post-thermostat NVE segment changed temperature by 0.084 K and total energy by 0.00038 eV. NVT still warmed by 8.51 K over its final 0.1 ps, and 8,667 mobile atoms moved over 1 Å relative to its start. Record this as an exploratory cold-baseline candidate, not a fully equilibrated static pass. The TTM grid smoke check passed with a combined energy residual of `1.48e-14 eV`. The finite-size LAMMPS/TTM link preflight stopped before any TTM step because the atom-to-grid phase volume exceeded the cell capacity (peak fill 1.316 at 10 Å); paired source/no-source coupling runs were not completed. See `interface_reconstruction_v6_exploratory/README.md` and `ttm/AGW_TTM_LINK_PREFLIGHT_V6.md`.
