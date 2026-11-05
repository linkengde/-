# P1G formal formal_case_A_marangoni_off_2100K_cut_cool_20261105

- Start from the common 300 K cold state and the same Ag70 QSGS geometry/input masks.
- Frozen mesh: 64x64x64; source power 8.931899 mW; Ag + Ti3SiC2 composite heating; negative d(sigma)/dT reference; solid drag 1e11 s^-1.
- Marangoni switch: off.
- Keep the normalized source fixed at 8.931899 mW. Stop heating when the measured Ti3SiC2 maximum reaches 2100 K, or when full configured input and bottom sink differ by no more than 0.01 mW for 10 consecutive steps. Do not throttle the source using a temperature predictor. Continue cooling and liquid Ag flow until the solidification criterion is met; retain 2300 K as the hard property ceiling.
- Full 3-D field snapshots every 0.05 ps. Full restart checkpoints every 0.2 ps, plus at heat-off and solidification.
- Safety horizon 1,000,000 ps is an emergency ceiling only; normal completion is the solver's solidification event. The shell watchdog is 72 h.
- No evaporation, recoil, ejection, production solver edits, or old Ag70 checkpoints.
