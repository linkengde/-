# Reproduction commands

Run from the P1H directory with the pinned Basilisk v2026-07-20 source set in `BASILISK_SRC` and the P1G input bundle available at `P1G_RUN_ROOT/inputs/basilisk`.

```bash
# Closed top, same 64^3 geometry and initial state
P1H_MAX_TIME_PS=0.16 \
P1H_MARANGONI=off \
P1H_ARC_PRESSURE_ENABLED=1 \
P1H_ARC_PRESSURE_PEAK_PA=1800000 \
P1H_ARC_PRESSURE_PROFILE=uniform \
P1H_ARC_PRESSURE_RAMP_PS=0.02 \
P1H_SOURCE_RAMP_PS=0.02 \
bash run_p1h_pressure_validation.sh smoke

# Open top; only the compiled top-boundary mode changes
P1H_MAX_TIME_PS=0.16 \
P1H_MARANGONI=off \
P1H_ARC_PRESSURE_ENABLED=1 \
P1H_ARC_PRESSURE_PEAK_PA=1800000 \
P1H_ARC_PRESSURE_PROFILE=uniform \
P1H_ARC_PRESSURE_RAMP_PS=0.02 \
P1H_SOURCE_RAMP_PS=0.02 \
bash run_p1h_pressure_validation.sh smoke-open-top
```

The open-top executable is built with `-DP1H_TOP_OPEN=1`. The 1.8 MPa value is only a software load for testing the code path, not an experimentally supported nanoscale pressure.
