# Model B: three-times source-flux sensitivity (0.25 ps)

Date: 2026-10-25. This is a source-amplitude sensitivity check on the exploratory rigid-obstacle Model B, using the paper-example peak flux of 15 eV/(ps·Å²), three times the project's nominal 5 eV/(ps·Å²). It is kept separate from the nominal-source pilot and is not a production parameter choice.

## Matched setup

The geometry, initial coordinates, 300 K velocity seed, timestep (0.5 fs), duration (500 steps), electron grid, Ag parameter tables, WCA obstacle, and fixed ceramic are identical to the nominal 0.25 ps pilot. Only the peak source amplitude changes. The Ag-masked source is not renormalized after applying the mask. The no-source control is the matched 500-step NVE baseline archived with the nominal pilot.

## Comparison

| Peak flux | Source power | Energy injected in 0.25 ps | Electron-to-lattice transfer | Source-weighted Ag lattice temperature at 0.25 ps | Maximum electron temperature | Combined-energy residual |
|---:|---:|---:|---:|---:|---:|---:|
| 5 eV/(ps·Å²) | 2934.740 eV/ps | 733.685 eV | 64.309 eV | 253.204 K | 1454.081 K | 0.058903 eV |
| 15 eV/(ps·Å²) | 8804.221 eV/ps | 2201.055 eV | 120.772 eV | 255.439 K | 2585.064 K | 0.059668 eV |

The three-times source raises the electron temperature substantially and transfers more energy to the lattice, but the reported source-weighted Ag lattice temperature changes by only 2.235 K at 0.25 ps. The maximum active-cell lattice temperature is 748.964 K; because it is a sparse-grid extreme, it is not a liquid-state criterion. No stable melt pool was demonstrated. The matched NVE drift is 0.057952 eV; the three-times-source residual differs from it by 0.001716 eV.

This outcome is consistent with the observed structural energy redistribution dominating early lattice-temperature response, but it does not isolate that mechanism from coarse-grained temperature noise or spatial heat diffusion. Do not infer from this test that Ag cannot melt at a longer duration or a different source. It does show that simply tripling the flux for 0.25 ps does not solve the present feasibility problem.

## Limits and reproduction

The B model still has no accepted cold baseline, physical Ag–Ti₃SiC₂ interface potential, mobile ceramic response, trajectory-based liquid classifier, or automatic heat-off controller. No cooling/solidification result was generated.

Reproduce the run from the repository root with the archived driver:

```bash
LD_LIBRARY_PATH=/workspace/.venvs/lammps/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} \
PYTHONPATH=/tmp/ag-ttm-mpich:/workspace/.venvs/lammps/lib/python3.12/site-packages \
OMP_NUM_THREADS=4 AG_TTM_PILOT_STEPS=500 AG_TTM_SOURCE_SCALE=3 python \
models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/geometry_constraint_B/b_ttm_reference_md_pilot.py
```

The captured stdout is `log.b_ttm_source_sensitivity_3x_0p25ps`.
