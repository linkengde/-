# P1H local arc-pressure validation

## Decision

The isolated pressure-load module and its basic VOF behavior are implemented and pass the checks listed below. This is **not** validation of a physical arc-pressure value and does not authorize the formal C/D heating runs. The existing gas-cap height has not been varied; open-versus-closed top-boundary behavior has now been compared over the initial-melt window. The height gate remains unresolved.

The work is on `codex/ag70-p1h-arc-pressure`, based on P1G commit `c262e2ea37f5b631da98a64b44604421150f36f4`. The P1G source and production enthalpy solver were not edited.

## Implemented load

The pressure is added to the existing CLSVOF interfacial normal-stress potential:

\[
\phi_{\mathrm{total}}=\sigma_{\mathrm{eff}}\kappa+p_{\mathrm{arc}}w_{\mathrm{Ag}}f_lw_{\mathrm{top}}.
\]

The mask is evaluated only in mixed VOF interface cells. Ag surface coverage is taken as `ag_marker/f`, the local Ag liquid fraction is applied, and `clamp(-n_z, 0, 1)` selects the upper-facing surface under the `f=1` sample / `f=0` gas convention. Positive pressure therefore gives a negative z force on the top surface. Pressure is not applied as a body force to every Ag cell or to the top domain boundary.

The default pressure switch is off. Runtime configuration supports a uniform pressure or a periodic-image Gaussian, an independent pressure center/radius, start/ramp/end times, and a peak pressure. A heat-source cutoff starts the same smooth off-ramp for the source and pressure; restart state records and checks the pressure/ramp configuration. Internal Ti3SiC2 pore surfaces are not represented by this external-surface VOF and receive no pressure load.

## Validation evidence

| Test | Result |
| --- | --- |
| Pressure helper/unit test | Passed profile, start ramp, mask, upper-surface sign, periodic Gaussian seam, and synchronized source/pressure ramp checks. Algebraic planar upper-face force was `-1.024e-13 N` for the test area and scale. |
| Uniform pressure on a planar VOF surface | Passed through dimensionless time `0.199900`. Center displacement was `2.42e-9` length units, volume drift `3.91e-9`, and final maximum speed `6.12e-10`; no resolvable deformation occurred. |
| Periodic-Gaussian VOF load | Passed in a synthetic nondimensional test with `sigma=1`, radius `0.18`, and `p=1.67` solver units (`Pi_p = pR/sigma ≈ 0.301`). At time `0.199900`, center-to-flank depression was `0.0242646` length units, greater than the half-cell threshold `0.015625`; volume drift was `7.56e-9`. This is only a force-path test, not an Ag pressure estimate. |
| Gaussian-load time-step sensitivity | At the same dimensionless end time `0.04`, maximum steps `1e-4`, `5e-5`, and `2.5e-5` gave center depression `0.00299097`, `0.00299565`, and `0.00299804`. The final two differ by `0.080%`; their measured maximum speeds were `0.146737` and `0.146548`. Maximum recorded relative volume drift was `2.33e-9`, `2.99e-9`, and `8.31e-10`, respectively. The response approaches a stable value as the step is reduced. |
| 64³ Ag70 geometry smoke run | Completed to the `0.16 ps` maximum-time stop with Marangoni disabled and a uniform `1.8 MPa` pressure used strictly as a software smoke-test input. The last fully logged heating state was `0.147782 ps`: maximum effective masked pressure `718,992 Pa`, integrated `F_z=-1.27423e-13 N` (downward), Ag-volume drift `0%`, peak skeleton temperature `1287.10 K`, liquid-Ag volume `0.190788 nm³`, and energy residual `3.73e-24 J`. A fixed-grid field snapshot was saved at `0.16 ps`. This run did not reach cooling or a peak-pool endpoint. |
| Pressure-configured restart | Restored the `0.16 ps` smoke checkpoint with matching pressure settings, recovered `Ag=4130.4257 nm³`, skeleton `1784.1705 nm³`, and source power `8.93190 mW`, then continued to the `0.161 ps` test stop. |
| Snapshot integrity | The HDF5 handoff is `64×64×64`, time `0.16 ps`, with 20 thermal/VOF fields and 4 pressure fields. Dataset values are finite; curvature uses Basilisk's finite `nodata` sentinel where no interface curvature is defined. |
| Closed/open top-boundary comparison | Completed matched 64³ runs to 0.16 ps with the same 11-cell (3.32 nm) gas cap, thermal source, initial state, and software-only uniform 1.8 MPa pressure. At the last logged liquid-Ag state (0.147782 ps), both runs had the same 0.190788 nm³ liquid volume, 1287.09664 K skeleton peak, 718992 Pa maximum effective pressure, -1.27423e-13 N integrated vertical pressure force, and 3.73e-24 J energy residual. At 0.16 ps the fixed-grid projected temperature, VOF, liquid fraction, material fields, and column-height map were bitwise identical. The maximum velocity difference was confined to the gas cap (3.12e-4 m/s); within the sample it was 2.92e-11 m/s, and within liquid/interface cells 2.33e-10 m/s. No resolvable surface depression was present. See `P1H_TOP_BOUNDARY_SENSITIVITY_REPORT.md`. |
| Cross-boundary restart guard | Passed: the open-top executable rejected the closed-top restart sidecar with exit code 2 and the expected boundary-mode mismatch message. |
| Open-top matching restart | Passed: the open-top 0.160 ps state resumed with matching settings and continued to 0.161 ps (18 steps, exit code 0). |

The synthetic Gaussian tests use nondimensional solver units. The `1.8 MPa` setting in the geometry smoke is not taken from a measurement or a macro-scale weld model and must not be described as a calibrated nanoscale arc pressure. The Ag surface-tension model remains an exploratory candidate, including its high-temperature endpoint treatment.

## Artifacts

- `p1h_arc_pressure.h`: pressure profile, mask, and time-ramp helpers.
- `p1h_free_surface_coupled.c`: isolated P1H solver copy; leaves P1G untouched.
- `p1h_pressure_vof_test.c`: uniform, local-load, and time-step VOF tests.
- `test_p1h_arc_pressure.c`: helper/sign/ramp unit test.
- `run_p1h_pressure_validation.sh`: reproducible build and validation entry point.
- `P1H_FIELD_OUTPUT_SCHEMA.md`: binary/HDF5 field and unit map.
- `convert_p1h_snapshot.py`: raw field-frame to compressed HDF5 converter.
- `results/dt_sensitivity_20261105/`: raw numerical logs and extracted comparison CSV.
- `results/geometry_preflight_0p16ps/`: compressed field snapshot, run histories, and logs.
- `results/restart_smoke_20261105/`: restart verification logs.
- `P1H_TOP_BOUNDARY_SENSITIVITY_REPORT.md` and `results/top_boundary_sensitivity_20261105/`: matched closed/open comparison, HDF5 projections, CSV, and run commands.
- `results/open_top_restart_smoke_20261105/`: same-mode restart log and compressed 0.161 ps checkpoint with sidecar.
- `results/restart_boundary_guard_20261105/`: cross-boundary restart rejection log.
- `results/validation_20261105T134109Z/`: pressure helper and VOF regression logs after adding the boundary mode.

Raw Basilisk restart dumps and uncompressed binary frames are retained in `/tmp/ag70-p1h/results/preflight_uniform_pressure_0p16ps_v3/` and `/tmp/p1h_raw_artifacts/`; interrupted diagnostic logs are preserved there as well. They are not committed because they are large. The committed HDF5 file is the postprocessed `0.16 ps` frame.

Build and test with the pinned Basilisk `v2026-07-20` source:

```bash
BASILISK_SRC=/path/to/basilisk/src bash run_p1h_pressure_validation.sh test
```

The input `1.8 MPa` geometry smoke can be reproduced with `P1H_ARC_PRESSURE_ENABLED=1`, `P1H_ARC_PRESSURE_PEAK_PA=1800000`, `P1H_ARC_PRESSURE_PROFILE=uniform`, `P1H_ARC_PRESSURE_RAMP_PS=0.02`, `P1H_SOURCE_RAMP_PS=0.02`, and `P1H_MARANGONI=off`; it is deliberately a short software smoke test.

## Required next checks before formal C/D runs

1. Vary the top gas-cap height while keeping the sample geometry, physical x/y period, and grid spacing fixed. The present solver is a cubic 64³ octree; increasing only z requires a rectangular-domain/mesh mapping or a carefully controlled supercell, so this isolated height comparison is not yet available. A 128³ supercell would also double the lateral periodic box and is not a one-variable test.
2. The closed-versus-ambient-pressure top comparison is complete for the existing 11-cell gas cap and initial-melt time window. It showed no measurable response within the sample, but does not establish cap-height convergence or behavior after substantial interface motion.
3. Quantify pressure-curvature force balance on the chosen interface and record grid/time-step behavior at the geometry scale.
4. Only then run pressure-on C/D from matched 300 K states. Keep pressure parameters labeled as sensitivity scenarios until a nanoscale-relevant source is established.

No formal C/D heating, Marangoni-plus-pressure comparison, cooling, or permanent-pit prediction was run in this phase. The current `0.302 nm` lattice is comparable to atomic spacing, while the flow calculation remains a continuum CLSVOF model; these outputs are exploratory.
