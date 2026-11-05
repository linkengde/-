# P1H top-boundary sensitivity check

## Scope

This was a short numerical boundary check on the existing Ag70 geometry. It compared the default closed/no-penetration top with an ambient-gauge-pressure outlet. Both calculations used the same 64³ geometry, 0.301857 nm grid spacing, 11-cell (3.320 nm) gas cap, initial state, thermal source, pressure input, and Marangoni-off setting. The pressure was uniform 1.8 MPa solely as a software smoke-test value; it is not a calibrated nanoscale arc pressure. Both runs stopped at 0.16 ps and did not cool.

The open-top mode sets zero-gradient velocity at the top and a zero-gauge pressure reference. The closed mode retains no-penetration velocity and the existing Neumann pressure condition. Snapshot case flags encode both Marangoni and top-boundary mode; restart state checks reject a checkpoint made with the other top-boundary setting.

## Results

At the last fully logged state with liquid Ag, 0.14778175 ps, the two runs reported identical values:

| Quantity | Closed top | Open top |
| --- | ---: | ---: |
| Peak skeleton temperature | 1287.09664 K | 1287.09664 K |
| Peak Ag temperature | 1244.79318 K | 1244.79318 K |
| Liquid Ag volume | 0.190788089899 nm³ | 0.190788089899 nm³ |
| Total Ag volume drift | 0% | 0% |
| Maximum effective masked pressure | 718992.109 Pa | 718992.109 Pa |
| Integrated vertical pressure force | -1.27422663403e-13 N | -1.27422663403e-13 N |
| Energy residual | 3.72736327974e-24 J | 3.72736327974e-24 J |

At the 0.16 ps projected fixed-grid snapshots, temperature, VOF, liquid fraction, material markers, curvature, enthalpy, and the derived column-height map were bitwise identical. The height map was computed as the z-integral of the projected VOF fraction multiplied by the grid spacing. Both maps remained flat at the output resolution, so no surface depression was resolved.

The velocity and pressure snapshots differ mainly in the top gas layer, where the boundary conditions act. Maximum velocity difference was 3.12e-4 m/s over the full box, but 2.92e-11 m/s in cells with sample VOF above 0.5 and 2.33e-10 m/s in liquid/interface cells. Maximum pressure difference was 14.09 Pa in the gas cap, 0.00782 Pa in sample cells, and 0.0103 Pa in mixed VOF cells. The histories through first melt likewise agree to the printed precision.

These results indicate no measurable closed/open top-boundary effect within the sample during the initial-melt window at the present gas-cap height. They do not establish boundary insensitivity after larger surface motion, later heating, or cooling, and they do not test gas-cap-height convergence.

The open-top checkpoint also resumed successfully with matching boundary and pressure settings from 0.160 ps to 0.161 ps (18 solver steps total, exit code 0). Conversely, attempting to use a closed-top checkpoint with the open-top binary was rejected before time integration.

## Gas-cap-height limitation

The current code and input use a cubic 64³ octree with equal cell spacing in all directions. The sample surface is at z=53 cells, leaving 11 gas cells to the top. Increasing only the gas-cap cell count while keeping the exact physical sample, x/y periodic length, and spacing requires a rectangular-domain or mesh-mapping change. A 128³ cubic supercell can add headroom, but it also doubles the lateral periodic lengths and raises memory, storage, and runtime; it would not be a one-variable test. That larger test was not launched.

The next valid height study should preserve the sample and x/y period while extending only z, then compare the same closed/open cases at the same physical time and pressure load. Until that geometry path exists, the gas-cap-height gate remains open and formal C/D runs remain unstarted.

## Reproducibility and files

- Solver code: `p1h_free_surface_coupled.c`; build open-top mode with `-DP1H_TOP_OPEN=1`.
- Converter: `convert_p1h_snapshot.py`; HDF5 attribute `top_boundary_open` records the mode.
- Numeric comparison: `results/top_boundary_sensitivity_20261105/top_boundary_comparison.csv`.
- Projected snapshots: `closed_top_0p16ps.h5` and `open_top_0p16ps.h5` in the same results directory.
- Open-top solver CSV histories: `results/open_top_pressure_smoke_20261105/`; projected field snapshots are in the top-boundary sensitivity directory.
- Same-mode restart smoke: `results/open_top_restart_smoke_20261105/`.
- Restart-boundary mismatch check: `results/restart_boundary_guard_20261105/`.

The HDF5 arrays are fixed-grid projections of the solver's adaptive fields; use the solver CSV for conservative liquid and total-Ag volumes. No physical crater, late-time pool, or solidification result is claimed here.
