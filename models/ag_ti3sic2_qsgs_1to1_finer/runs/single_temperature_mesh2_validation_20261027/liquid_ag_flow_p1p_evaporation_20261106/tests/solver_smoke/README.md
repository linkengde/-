# P1P integrated source smoke test

These files are software diagnostics, not Ag evaporation results. The smoke
run used the existing 64 x 64 x 64 Ag70 material map plus a 12-cell open gas
cap, but it used `synthetic_property_table.csv` and a manufactured 1400 K
initial thermal field. Its synthetic saturation pressure and latent heat must
not be interpreted as Ag properties.

The `P1P_TEST_STOP_AFTER_SOURCE` build ran the actual top-connected-gas search,
PLIC area reconstruction, interface-temperature interpolation, Hertz-Knudsen
source, and conservative cell update once, then exited before the inherited
thermal solve. At `t = 0.00181818182 ps`, it found `2.98320138895e-16 m^2` of
exposed liquid-Ag PLIC area, a maximum synthetic flux of `1.55501282076
kg/(m^2 s)`, and removed `1.25317752386e-30 kg`. The reported Ag mass residual
was `-4.5838895855e-13%`. The exact run log and source CSV are included here.

A second run without the stop-after-source test hook proceeded into the
thermal solve. The P1P source ran, but the inherited nonlinear enthalpy solve
did not converge from the manufactured hot initial condition: at the same
time it stopped after 81 iterations with a `48.9522415 K` update. That failure
is retained in `synthetic_thermal_failure.stderr`; it means the complete
thermal-flow runtime path has not passed. No restart, long run, or physical
Ag70 evaporation calculation was performed.

With a 1 fs maximum-step setting, a further P1P run completed one full
source-plus-thermal step. Its Ag-volume drift was `-1.84719291e-10%` and its
global energy residual was `-5.97478060099e-24 J`. The next step failed in the
thermal nonlinear iteration at `0.00026446281 ps` after 81 iterations with a
`102.081935 K` update. A matched run with evaporation disabled failed at the
same time, iteration count, and update, showing that this diagnostic failure
comes from the inherited thermal solve under the manufactured hot initial
condition, not from P1P. Both runs' logs and CSVs are retained.

Both runs used the Basilisk `v2026-07-20` source snapshot recorded in the
project manifest and the same P1M input directory. Rebuilds must set
`BASILISK_SRC` and compile with `P1H_RECT_NZ=76` and `P1H_TOP_OPEN=1`. The
stop-after-source macro is test-only and must not be used for production runs.
