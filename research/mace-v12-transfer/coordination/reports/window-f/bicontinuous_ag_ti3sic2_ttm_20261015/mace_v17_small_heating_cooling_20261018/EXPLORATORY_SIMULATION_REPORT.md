# Ag–Ti₃SiC₂ MACE V17 exploratory heat/cool run

## Result

The requested CPU-limited exploratory run was carried through local heating and a separate heat-off/cooldown segment. It produced an OVITO-ready 30-frame trajectory through 0.13 ps, but **did not establish an Ag melt pool or resolidification**. At the highest sampled local Ag temperature, the 92-atom hot region remained predominantly FCC-like. After heat-off, its local temperature fluctuated strongly over only 0.02 ps; that interval is too short to demonstrate complete cooling.

The result is useful as an early trajectory and a throughput/stability screen. It is not a melting or welding result.

## Model used

- The original 39,486-atom Scheme A model and its V17 frozen static force result were checked by SHA before model work. Source data SHA256: `e2b62a68c46e23037370a6ca20951aa0774edea2d728d3c3aadbcba4b3755674`.
- Because that full structure's V17 run-0 already used about 29.04 GiB RSS in a 32 GiB container, full-system dynamics were not attempted. The available environment has 4 CPU-equivalents, 32 GiB RAM, and no CUDA device.
- An independent 1,926-atom layered interface slab was built from ideal COD 9009647 Ti₃SiC₂ and ideal FCC Ag(111), not from any melt-pool trajectory. It contains 1,200 Ti₃SiC₂ atoms and 726 Ag atoms, or **66.6749 wt% Ag**. This differs from the original approximately 59.95 wt% composition and from the original bicontinuous topology.
- The small slab has a C-terminated Ti₃SiC₂(0001)/Ag(111) interface, XY periodicity, free Z surfaces and 18 Å top vacuum. The original full model has free X and periodic Y/Z, so the slab is a representative interface model rather than a boundary-condition or morphology copy. Only the bottom 100-atom Ti plane is fixed; the remaining Ti₃SiC₂ backbone can vibrate and conduct heat.
- Ag was compressed in-plane by 3.8907% to match the finite cell. Initial cross-phase minimum separation was 2.7048 Å. No atoms were deleted to suppress contacts.
- Only the frozen MACE V17 weight was used: SHA256 `559b12c024fe8e2155e26ffc7ca48cf9fa1c2682ea37ff9c5367af14a5cc43d3`; CPU float64, ASE/MACE 0.3.16, Torch 2.5.1+cpu. No EAM, TSCBOP, or Morse terms were added.

## Calculation record

- Initial MACE static evaluation was finite. A 60-step FIRE screening relaxation reduced Fmax from 2.731 to 0.897 eV/Å. A partial 300 K baseline ran 0.02 ps with no NaN, atom escape, or numerical collapse; it was not a full stability validation.
- Local heating used a 1,900 K Langevin source on the top three Ag(111) planes within an 8 Å lateral radius of the cell center. The original hot set was 92 Ag atoms. Mobile Ti₃SiC₂ atoms in the support-side sink layer were coupled to 300 K; other mobile atoms evolved without a global thermostat.
- The initial 1 fs heating segment reached 0.09 ps. At its last frame Fmax was 6.384 eV/Å, with a 2.243 Å Ag–Ag contact and a 1.969 Å Ag–TSC contact. Heating stopped at that point as a precaution. A separate continuation started from the saved positions and velocities at 0.2 fs and ran 100 steps to 0.11 ps. Fmax fell to 4.037 eV/Å; the corresponding minimum Ag–Ag, Ag–TSC, and Ti–C distances were 2.295, 2.013, and 1.832 Å. The high force did not disappear immediately at the restart, but it did not run away. This was not a paired same-state 1 fs/0.2 fs comparison, so the change cannot be attributed solely to timestep.
- The continuation script reselected the circular hot region at the 0.09 ps restart, resulting in 90 thermostatted Ag atoms rather than the original 92. The final offline trajectory analysis uses the fixed original 92 atom IDs for consistent temperature, order, and mobility comparisons. The cooling diagnostic script also reselected its spatial window after atoms had moved; this affected only its raw diagnostic atom count, not the dynamics or heat-removal bath. The published OVITO metadata recomputes temperatures on the same original 92 IDs.
- A separate 0.2 fs cooling segment then ran 100 steps to 0.13 ps with the Ag source off and the 300 K Ti₃SiC₂-side heat sink active. The 92-atom fixed-ID hot-region temperature was noisy and nonmonotonic over this 0.02 ps interval. This is a cooling-start trajectory, not a completed quench.

## Structure and phase evidence

Ag order was measured with per-atom Ag–Ag Steinhardt q₆ using a 3.45 Å cutoff, alongside coordination, displacement, and connectivity. q₆ below 0.35 was used only as a local disorder tag; a liquid candidate also required a large connected low-q₆ cluster and substantial Ag migration.

- Baseline hot-region median q₆: 0.5724.
- Hottest sampled frame: 0.08 ps, 1,624 K; hot-region median q₆ 0.5286, median coordination 12, zero of 92 hot Ag atoms with q₆ < 0.35, and RMS displacement 0.275 Å from the 300 K reference.
- Final cooldown frame: 0.13 ps; the same 92 Ag atoms had median q₆ 0.5034, median coordination 12, one atom below q₆ = 0.35, largest connected low-q₆ cluster of one atom, and RMS displacement 0.361 Å.
- The multi-indicator melt criteria were not met. **Melt-like connected Ag region: 0 atoms; size and dimensions: none identified.** The one low-q₆ Ag atom in the final frame is an isolated local fluctuation, not a liquid region.
- The Ti₃SiC₂ neighbor graph remained one component at a 3.25 Å cutoff in every analyzed frame. Minimum backbone connectivity fraction was 1.0. At the hottest frame mobile-TSC RMS displacement was 0.259 Å; the minimum sampled Ti–C distance was 1.814 Å. This indicates a connected vibrating backbone over the short trajectory, not validated long-term stability.
- Since melting was not established before heat-off, **resolidification evidence is also not established**. No reliable weld-strength or joint-formation conclusion follows from these data.

## Cost and limitations

The three thermal segments (0.09 ps at 1 fs, followed by 0.02 ps at 0.2 fs heating and 0.02 ps at 0.2 fs cooling) used about **5.33 aggregate CPU-hours** and **88.3 minutes wall time**, with peak RSS about 17.7 GiB. The first heating segment used about 28.5 minutes wall time; each 0.02 ps 0.2 fs segment took about 30 minutes. On this representative slab, a linear extrapolation from the 1 fs heating segment is about 5.3 wall-hours per 1 ps; it is a rough budget estimate, not a performance guarantee. At 0.2 fs, the measured pace is about 25 wall-hours per physical ps.

LAMMPS-MACE acceleration was investigated but not accepted: the CPU ML-IAP path lacked the MACE ghost-feature exchange callback, and the separately built Kokkos path aborted with an invalid-pointer error before producing forces. The original installed LAMMPS binary lacks ML-IAP. ASE/MACE was therefore used for the actual trajectory. OVITO is not installed in this container; the final file is standard Extended XYZ and was not locally rendered.

The 39,486-atom static V17 calculation's ~29 GiB RSS leaves inadequate headroom for full-size dynamics in this 32 GiB environment. No full-size MD, DFT, TTM, or MACE training was started. Original structure, weights, and archived runs remain unchanged.

## Final flags

```text
MACE_V17_SMALL_MODEL_READY=YES
MODEL_ATOM_COUNT=1926
300K_BASELINE_STATUS=PARTIAL_0.02PS_NO_NUMERICAL_COLLAPSE_NOT_VALIDATED
HEATING_STARTED=YES
COOLING_STARTED=YES
AG_MELTING_EVIDENCE=NOT_ESTABLISHED
AG_RESOLIDIFICATION_EVIDENCE=NOT_ESTABLISHED
OVITO_TRAJECTORY_PATH=AG_TSC_MACE_V17_HEAT_COOL_OVITO.extxyz
ACTUAL_CPU_TIME=~5.33 CPU-hours for the three heat/cool segments
FULL_39486_ATOM_MD_FEASIBLE_IN_CURRENT_ENVIRONMENT=NO
```
