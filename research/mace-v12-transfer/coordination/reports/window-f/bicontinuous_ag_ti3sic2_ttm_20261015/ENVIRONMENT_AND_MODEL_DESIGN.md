# Ag–Ti₃SiC₂ bicontinuous TTM melt-pool model: file audit and design

**Status: preparation only.** This report records the current F cloud environment, files actually present or recoverable from GitHub, and a proposed small-to-large model route. No atomic structure was generated, no potential was changed, and no LAMMPS/DFT/TTM calculation was started.

## Environment and repository access

- Working tree: `/workspace/-`, repository `https://github.com/linkengde/-.git`.
- This task is isolated on branch `cloud-f-bicontinuous-ttm-pool-20261015`, based on `origin/main` at `e91ef26d175072553b409e08f5d9a826929078c7`. The prior `cloud-f-tscbop-refit` branch was left intact.
- GitHub fetch succeeded for the Stage67 handoff and model-build references listed below. Their files were not assumed to be present in this checkout.
- Current resources: 4 CPU-equivalents by cgroup quota (`cpu.max=400000 100000`), 5 logical CPUs exposed (`0–4`), 32 GiB memory limit, about 27 GiB free disk and about 2.0 million free inodes.
- `/usr/bin/mpiexec` is present. No `lmp` executable or LAMMPS Python module was found in the current environment, so even the requested static `run 0` is unavailable here. No packages were installed.

## Potential files

These three files are present in the current checkout. Their bytes and hashes match the Stage67 handoff inventory. No replacement was selected.

| File | Size | SHA-256 | Identity / status |
|---|---:|---|---|
| `research/mace-v12-transfer/coordination/reports/window-b/current_heating_potential_handoff_20261009/potentials/Ag_u3.eam` | 36,854 B | `cb3d605cd27ec6efb50c4dd7462c2ab3d967e63775dc4400de9134315a44d807` | Foiles universal-3 Ag EAM; FCC lattice parameter 4.090 Å in file header. |
| `research/mace-v12-transfer/coordination/reports/window-b/current_heating_potential_handoff_20261009/potentials/TSCBOP.tersoff` | 3,330 B | `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7` | Plummer–Tucker Ti–Si–C bond-order file. The archived fixed-bulk comparison reports about 1.99 eV/Å force RMSE; its report does not qualify it as a production skeleton potential. |
| `research/mace-v12-transfer/coordination/reports/window-b/current_heating_potential_handoff_20261009/potential.interface_calibration_v0.inc` | 641 B | `34deced5e08174ba369290b87f5c3bfd2b284b272991adf302b2b0d28c12815f` | Ag–Ti/Si/C Morse calibration-v0 include; explicitly not production-validated. |

The cross terms in calibration-v0 are Ag–Ti `(D0=0.94814194 eV, α=2.0 Å⁻¹, r0=2.33708864 Å)`, Ag–Si `(0.93944948, 2.0, 2.40407837 Å)`, and Ag–C `(1.35077230, 2.0, 2.15255924 Å)`, each with 8 Å cutoff. The accompanying README says they were fitted to a single-plane HCP surrogate and warns that the old heating branch was paused over cross-potential energy release. These are preserved only as the historical candidate identity; the file is not evidence of reliable 3D interface contacts or solidification behavior.

The fixed-bulk force warning is documented in E commit `7701ecfdd4e565c4d4be16b889d091c44aa8b648`, file `research/mace-v12-transfer/coordination/reports/window-e/low_mismatch_217_dft_lite_20261013/minimal_135atom_model/potential_development_decision_20261014/TSCBOP_force_error_environment_analysis.md`: the same TSCBOP gives 1.9965/1.9872 eV/Å total force RMSE on one fixed 48-atom parent structure at two k meshes. This is a warning about the current skeleton potential, not a new calculation in this task.

## Stage67 and Ce file availability

The current checkout does not contain the Stage67 TTM data tree. It is recoverable from GitHub branch `cloud-d-stage67-300x300-extension-handoff-20261013` (commit `f4faca36e3987ddc577c49c6ce8c27a1731073da`), under `stage67_300x300_extension_handoff_20261013/archive_parts/`. I reconstructed its 41 ZIP parts in `/tmp`, verified the resulting ZIP SHA-256 `6cbef8dbf96ee2f94bf5c2929b0b622902db12898f879cde00e6c97c56a59547`, passed `unzip -t`, and verified all 48 entries in the package SHA-256 manifest. This was a temporary read-only copy; the old project tree was not altered.

| Stage67 asset in the verified package | SHA-256 | Finding |
|---|---|---|
| `input/implementation_reference/Stage67_200A_1fs_driver.py.reference.txt` | `633a8a6a4d4927df8ee307f12febc4d8b8c18f6430b062d54fb43cd0507dbb2b` | 78,138-byte implementation reference with absolute Windows paths; not a standalone runnable input. Its older driver records `SOURCE_FACTOR=4` and 64,000 eV/ps, so it must not be copied into the requested new pulse. |
| `input/ttm/material_tables/stage43_Ag_Ce_table.csv` | `71591dc3f41641ca1ddc0e1b6fabf9c363614a6659e267e74e5fd8db09ec4679` | Ag Ce table, 300–20,000 K, source column says Medvedev author dataset. |
| `input/ttm/material_tables/TSC_Ce_Ue_hybrid_300_20000K.csv` | `750519ecb9d62870cbbf69e47a399039ccd1b63f954d37cc4c396582f23c7328` | Ti₃SiC₂ Ce/Ue table, 300–20,000 K; its upper-temperature rows are tagged `DOS_EXTRAPOLATION`. |
| `input/ttm/material_tables/G_TSC_BASE.csv` | `9de12c26f3239fc277e9911f8cd702d362e9e0ddd0f25f6775c044ac53250d01` | TSC coupling table. |
| `input/ttm/material_tables/stage43_Ag_kappa_table.csv` | `dbbe4b7d9ae1d63a2728b5b728c99f56ed8953b80e98b4b190b2b17eee632ce` | Ag thermal conductivity table. |
| `input/ttm/material_tables/TSC_kappa_BASE.csv` | `5f29d1da5acecbe646619348748db61b42701904fac29cb3f60f49179dcc302` | TSC thermal conductivity table. |
| `input/ttm/200A_spatial_maps_NOT_FOR_300A/stage67_source_energy_audit_4x.csv` | `d9e52428e4677e6fed475634f7124a4f93593b35794d77893b23801b2ebee4420` | Old 200 Å mapping only; not reusable as the new map. |

The accompanying old-grid heat-source audit describes a finite-duration TTM source entering the electron subsystem first, with an exponential depth weight and circular channel. On its 10 Å grid it used `R_arc=30 Å`, `delta=20 Å`, `I0=5 eV/(ps·Å²)`, discrete arc area 3,200 Å², and 16,000 eV/ps single-source power. The old 4× driver multiplied this to 64,000 eV/ps. The discrete normalization must be rebuilt and checked on the new mesh; the intensity is an initial test value, not a measured arc power. The new request does not yet specify pulse duration, source center, or whether the 150 Å X length includes vacuum.

The package reports these missing Stage67 inputs: the manifest-referenced `in.regions_x_surface.inc`, a standalone LAMMPS input matching the new geometry, a target-specific 15×30×30 spatial map, and the exact Stage63 seed `stage63_equilibrium_100.data` (expected hash `e193a7cde41a06fa04fd32ea9a80c60c4d5089d2a036a0c7fe6957ca57b4813f`). That Stage63 seed is not needed for the newly proposed crystal/TPMS build, but the missing region include and standalone implementation prevent direct Stage67 input reuse. The archived step-1650 structure is a 150.1×200×200 Å fully periodic evolved melt-pool state; it is not the requested new crystalline model and will not be used as its geometry.

## Proposed bicontinuous model

Use a periodic gyroid level-set as the phase selector while building each phase from its own perfect crystal lattice. Let

`g(r) = sin(kx)cos(ky) + sin(ky)cos(kz) + sin(kz)cos(kx)`, where `ki = 2πri/p`.

Assign one side of the surface to Ag and the other to Ti₃SiC₂; tune the threshold to the target Ag fraction, then tune it against the actual retained atom masses. A starting morphology period of about 50 Å gives six periods across each 300 Å periodic direction. The target phase volumes should not be set to 60/40: using the Ag EAM lattice density (~10.47 g/cm³) and the archived PBE Ti₃SiC₂ cell (~4.56 g/cm³), a 60/40 mass ratio corresponds to about **39.5 vol% Ag and 60.5 vol% Ti₃SiC₂**. This volume split is only an initialization estimate; the final mass ratio comes from atom counts.

For low-strain Y/Z periodic matching, a useful orientation candidate is Ti₃SiC₂ `c` along X, with its perpendicular basal supercell vectors of lengths `a=3.0575 Å` and `√3a=5.2957 Å` along Y/Z. Orient FCC Ag with `[001]` along X and `[110]` / `[1̅10]` along Y/Z, whose translation lengths are `a_Ag/√2=2.8921 Å`. Near the 300 Å targets, 104 Ag repeats span 300.775 Å; 98 and 57 Ti₃SiC₂ repeats span 299.635 Å and 301.857 Å. Fitting those crystals into exact 300 Å box vectors takes less than 0.62% in-plane strain per phase. A common-cell strain and its residual stress still need explicit review before use.

For the 20k–50k representative, start near a **60 Å solid thickness × 98.2 Å × 101.1 Å** material region, about **46,000 atoms** at the estimated 60/40 mass ratio. If using the Stage67 surface convention, add a 50 Å X vacuum margin before the material, making the box X length about 110 Å while retaining roughly 46,000 atoms. This small sample has two approximate morphology periods across Y/Z and is for local geometry/force/thermal stability checks only. Its 30 Å-radius liquid zone would have only about 40 Å clearance between periodic images, so it cannot certify the no-periodic-spanning condition for the full model.

For the full target, first freeze the X interpretation. The Stage67 convention places the free surface at X=50 Å and material on the +X side. If the requested 150 Å is the whole box, that leaves 100 Å of material and the estimate is about **697,000 atoms**. If 150 Å means material thickness in addition to 50 Å vacuum, the box must instead be about 200 Å in X and the estimate is about **1,045,000 atoms**. These estimates use unchanged crystal number densities and the same 60/40 mass split; they are not generated atom counts. The exact target-Y/Z cell with a 30 Å-radius pool gives 120 Å from the pool edge to a face and 240 Å between repeated pool edges, but actual liquid-cluster percolation must still be measured.

Construction and acceptance order:

1. Create full FCC Ag and Ti₃SiC₂ crystal lattices in a common cell with `pbc=(false,true,true)`, with separate crystal origin shifts. Apply the gyroid mask, then adjust the mask threshold to match retained atom masses.
2. Search deterministic origin shifts for the smallest Ag–Ti/Si/C contacts. Apply no random displacements. Reject any atom-level culling that breaks either phase's connected graph; do not invent a contact cutoff from the unvalidated Morse file.
3. Check per-phase largest connected components and spans in X/Y/Z (with wrap only in Y/Z), formula counts, actual mass fraction, cell and PBC, all pairwise minimum distances, and the morphology/atom mapping hashes. Require the solid Ag network and Ti₃SiC₂ skeleton each to remain connected; separately track the liquid-like Ag cluster after heating and reject periodic Y/Z wrapping.
4. Only after obtaining a LAMMPS build and reviewing the force-field applicability, run the small model in order: static `run 0`, energy minimization, 300 K equilibration, then source-free NVE drift test. Log maximum force, energy, element- and phase-resolved temperature, nearest contacts, and phase connectivity. If temperature rises, isolate Ag-only, TSC-only, and cross-interface terms on copied diagnostic inputs before changing any parameters.
5. After small-model checks pass, create the full cell and new TTM maps from its actual coordinates. For a 150×300×300 Å box and 10 Å cells the nominal mesh is 15×30×30 (13,500 cells); rebuild material fraction, `G_mix`, initial Te, and source mask. Recheck discrete source power, pulse duration, source center, surface orientation, and energy closure before any short heating test. Do not proceed directly to 10 ps or 100 ps.

## Open items before any structure generation or calculation

1. Confirm whether X=150 Å includes the suggested 50 Å vacuum or is material thickness.
2. Freeze pulse center, pulse duration, the heated X face, and nonperiodic boundary behavior; the available old driver cannot decide these for the new experiment.
3. Review the TSCBOP bulk-force mismatch and independently validate the Ag–TSC cross interaction. Current files are available, but not certified for the target melt/freeze process.
4. Provide/install an approved LAMMPS executable for this environment before static or dynamic validation; current cloud has no `lmp`.
5. Audit the Ce/Ue, kappa, and G tables' units, interpolation, and valid temperature range. Avoid relying on the tagged `DOS_EXTRAPOLATION` region without a stated validity decision.

The present deliverable is only an availability audit and geometry design. No model file, trajectory, restart, or TTM map has been created.
