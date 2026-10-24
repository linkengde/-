# Fixed-scaffold Ag cold pilot

This is a **no-source** 300 K early-stability check, not a melt-pool calculation or thermal-equilibrium test. All Ti/Si/C atoms and the original bottom layer are fixed; mobile Ag atoms are thermostatted at 300 K with a 0.1 fs step for 0.005 ps (50 steps). It uses the repository's current EAM + TSCBOP/Tersoff + calibration-v0 Morse interactions.

The initial coordinates come from the existing Stage 1 fixed-skeleton Ag-only FIRE minimization (`testA_ag_only_fire`, 200 iterations, `dmax=0.01 Å`). The compressed snapshot contains only the 180,363 mobile Ag atoms. Its uncompressed SHA-256 is `8bb02ecb7f12b493e3921df211c7038540098ffb6e239e9692bd0da07a137e2d`; the archived gzip SHA-256 is `b4254cd8c7b196a9588bf3b3ae4729743d0788fb7e8c83c1e70e435a946e6333`.

The input defines the original fixed/bottom/mobile groups before reading the relaxed Ag coordinates. This preserves layer membership even though a few relaxed Ag atoms may cross the original z=20 Å group boundary. The endpoint CNA calculation runs only after the MD segment because CNA is unusually expensive on this high-defect, high-coordination initial geometry.

Run from the model directory:

```bash
lmp -in cold_pilot/in.fixed_skeleton_ag_cold_0p005ps_from_step200 -log cold_pilot/log.fixed_skeleton_ag_cold_0p005ps
```

The result is conditional on the frozen-scaffold approximation and current potentials. Fifty steps only screen immediate instability; they do not validate a weakly mobile backbone, thermal equilibrium, or authorize local heating. Assess Ag temperature/force drift and the backbone's physical forces before using it as a starting state.

The 50-step trajectory and released-scaffold `run 0` logs are archived as `log.cold_pilot_0p005ps` and `log.physical_force_run0`. Compressed mobile-Ag coordinate snapshots at the initial relaxed state and the 50-step endpoint are archived alongside them.

For an early weak-restraint comparison from the same Ag-relaxed snapshot, run the parameterized `in.spring_cold_0p005ps` with `-var K 1` and `-var K 10`. These 50-step runs keep the bottom fixed, integrate the mobile region with NVE, thermostat the original bath layer, and restrain the mobile Ti/Si/C atoms with `spring/self`. At each endpoint the spring and dynamical fixes are removed before a `run 0` force check. They are short transient screens only; they do not establish equilibrium.

The K=1/K=10 comparison, result logs, and compressed endpoint Ag snapshots are archived in [RESULTS_WEAK_SPRINGS_20261024.md](RESULTS_WEAK_SPRINGS_20261024.md).
