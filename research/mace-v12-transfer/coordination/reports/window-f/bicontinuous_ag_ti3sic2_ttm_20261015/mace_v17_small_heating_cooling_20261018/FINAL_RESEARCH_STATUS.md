# Final status

```text
MACE_V17_SMALL_MODEL_READY=YES
MODEL_ATOM_COUNT=1926
300K_BASELINE_STATUS=PARTIAL_0.02PS_NO_NUMERICAL_COLLAPSE_NOT_VALIDATED
HEATING_STARTED=YES
COOLING_STARTED=YES
AG_MELTING_EVIDENCE=NOT_ESTABLISHED
AG_RESOLIDIFICATION_EVIDENCE=NOT_ESTABLISHED
OVITO_TRAJECTORY_PATH=AG_TSC_MACE_V17_HEAT_COOL_OVITO.extxyz
ACTUAL_CPU_TIME=~5.33 CPU-hours across heat and cooldown stages
ACTUAL_WALL_TIME=~88.3 minutes across heat and cooldown stages
FULL_39486_ATOM_MD_FEASIBLE_IN_CURRENT_ENVIRONMENT=NO
```

The delivered trajectory is 30 frames of 1,926 atoms, from 0 to 0.13 ps, and includes a 0.02 ps heat-off segment. The hot Ag region remained mostly ordered: at the hottest sampled frame its q₆ median was 0.529, coordination median 12, and 0/92 atoms were tagged q₆ < 0.35. The final cooling frame had 1/92 atoms below that tag, with no connected multi-atom low-order cluster. No melt region was identified; no resolidification claim is supported.

The Ti₃SiC₂ backbone stayed fully connected at the 3.25 Å graph cutoff throughout the saved frames. The local cooling trace was noisy and nonmonotonic at 0.02 ps, so it is a cooling start rather than a demonstrated completed cooldown.

Full-size 39,486-atom dynamics were not run: its earlier V17 frozen static inference already used about 29 GiB RSS under a 32 GiB limit, with no GPU available. The 1,926-atom slab differs from the original in composition, topology, and boundary orientation; details and all limitations are in `EXPLORATORY_SIMULATION_REPORT.md` and `MODEL_DESIGN_AND_COST.md`.
