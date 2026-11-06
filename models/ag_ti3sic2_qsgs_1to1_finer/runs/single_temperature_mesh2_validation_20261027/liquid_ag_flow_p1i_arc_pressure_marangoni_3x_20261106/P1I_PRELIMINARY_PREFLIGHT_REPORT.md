# P1I preliminary preflight (0.16 ps)

This archived run is an initialization and short-time coupled-path check only. It did **not** reach the 2100 K hold trigger and does not establish the outcome of the 3 ps hold, cooling, or solidification run.

## Configuration

- Grid: 64 × 64 × 76; the Ag70 geometry and 53-cell sample region are unchanged, with a 23-cell gas cap.
- Heating: existing composite source, 8.93189911185 mW reference power, 0.02 ps ramp.
- Marangoni: enabled; multiplier 3 applies only to the existing tangential thermal acceleration. Normal capillary stress is unchanged.
- Arc pressure: periodic Gaussian, configured peak 501 MPa, 5 nm radius, 0.02 ps ramp.
- Ti₃SiC₂ property table: 300–2300 K; 2100 K hold trigger and 2300 K hard property ceiling.
- Time step cap: 20 fs; five OpenMP threads; stopped at 0.16 ps by the preflight limit.

## Observed results

- Maximum skeleton temperature at stop: 1287.09 K. The 2100 K trigger was not reached.
- First recorded nonzero liquid Ag volume: 0.01068 nm³ at 0.13556 ps; at 0.14778 ps it was 0.19077 nm³.
- Maximum recorded speed at 0.14778 ps: 2.0048 m/s.
- Configured pressure was 501 MPa; maximum effective interface pressure recorded by 0.14778 ps was 184.77 MPa, with integrated z force −3.32167 × 10⁻¹¹ N.
- Ag volume drift stayed at floating-point noise (about −1.78 × 10⁻¹⁴%).
- At 0.14778 ps, the reported energy residual was −2.33 × 10⁻²² J.
- The raw pressure-projection residual reached 2.677 × 10²² after liquid flow began. This is recorded as a diagnostic quantity, not claimed as a converged error norm; compare it against the prior 167 MPa path and inspect solver behavior in the full run.

The exact run logs, CSV histories, compressed initial and 0.16 ps field snapshots, configuration, and SHA-256 manifest are in `results/preflight_0p16ps_20261106/`.

## Next run and stopping logic

The full runner is `run_p1i_3x_hold3ps.sh`. Before the 2100 K hold begins, it now permits heat-off if input power and bottom-sink power differ by at most 0.01 mW for 10 consecutive steps. If the skeleton reaches 2100 K first, heating and pressure remain active for the full 3 ps hold unless the predictive 2300 K property guard or divergence guard intervenes. Heat and pressure then ramp off together; the run cools until the solver reports solidification. The runner retains the two latest valid restart checkpoints and compresses snapshots while the run proceeds.

RAM was not the limiting resource at preflight: approximately 15 GiB was available. The nearly full `/tmp` and `/workspace` filesystems were the constraint, so redundant older checkpoints from completed runs were pruned while preserving each run's final state and latest valid restart.
