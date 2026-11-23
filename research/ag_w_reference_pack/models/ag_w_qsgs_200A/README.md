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

TTM heating remains the intended heating method. The stock `fix ttm/mod` does not support this free-surface boundary. The legacy Ag-only spatial source has been mapped onto this model's actual Ag atoms in `ttm_source_preflight/`; at `S(t)=1`, its Newton-layer Ag power is 9,784.984 eV/ps. The mask uses local atom positions, not the global Ag volume fraction. Pulse timing, W electronic-property data, and a validated mixed-phase closure are still missing. No thermostat pulse is substituted for TTM. See `ttm_source_preflight/TTM_SOURCE_PREFLIGHT.md` and the project-level `../../TTM_HEATING_READINESS.md`.
