# P1P follow-up: distinguish the evaporation-off baseline failure from C-case stop

Date: 2026-11-07
Branch: `codex/ag70-p1p-evaporation`
Reviewed baseline commit: `c256cfe2ba4345027a1ce6c4e4ae171452047f0c`

## Decision

Formal flow/heating remains frozen. No production solver, geometry, heat source, pressure, or Marangoni settings were changed, and no new run was launched for this audit.

The previously documented second-step failure in the evaporation-off synthetic smoke test is real and must be treated as an inherited baseline solver failure. It is independent of evaporation. The later C-case stop is a different immediate failure mode: an evaporation property-range rejection in a few nearly empty interface cells, before the thermal nonlinear solve begins. Both paths use the advected enthalpy/composition state and therefore justify a shared transport/inversion audit, but the available evidence does not prove one common root cause.

## Evidence and failure locations

| Record | State and setup | First successful step | Stop | What it establishes |
|---|---|---|---|---|
| `small_dt_no_evap_control` | Synthetic software fixture; 64×64×76; manufactured 1400 K initial thermal field; synthetic properties; 1 fs maximum step; evaporation off | `t=9.09090909e-05 ps`; `Tmax=1409.79383 K`; 16 nonlinear iterations; final update `0.00502619659 K` | At `t=0.00026446281 ps`, thermal nonlinear solve did not converge after 81 reported iterations; update `102.081935 K` | The inherited thermal solve can fail without evaporation. This is a software diagnostic, not a physical Ag result. |
| `small_dt_evaporation` | Same synthetic setup and initial state; evaporation on with a synthetic two-row property table | Same time, peak temperature, nonlinear iteration count and update as the no-evaporation control to displayed precision | Same thermal-solver stop time, reported iteration count and update | Enabling the source did not cause or prevent this observed second-step failure. The first-step enthalpy difference was about `2×10⁻²⁸ J`; the energy residuals differed by about `3.7×10⁻²⁹ J`. |
| Original P1M no-evap control | Same synthetic smoke family, no evaporation | Same first-step state and convergence record | Same second-step stop | Independently preserves the inherited no-evap failure in the original control record. |
| Full C case | Exploratory physical-property run; evaporation enabled | Continued through `9.30330989 ps` | Seven cells rejected by P1P property/composition validation before the thermal update; no cooling followed | This is not the second-step nonlinear-solver event. C remains incomplete and has no solidified crater assessment. |
| Instrumented C replay | Restarted at `9.00332084 ps` from the C state | Replayed to `9.28035319592 ps` | Five cells returned `P1P_ERR_RANGE`; reconstructed interface temperatures `4925–5320 K` versus table range `1234–2436 K` | Localizes a similar immediate property-range failure, but it stopped `0.023 ps` earlier and found a different cell count; it is not an exact reproduction or an upstream root-cause diagnosis. |

The source-only smoke intentionally exits after the P1P interface-source update and before the inherited thermal solve. It verifies that the source path executes; it does not validate the complete coupled runtime.

## Shared code path and current uncertainty

The solver initializes a bulk volumetric enthalpy field `Hvol`. Basilisk advects `Hvol` as a general tracer, while Ag, Ti₃SiC₂, and pore fractions are attached to the composite VOF field `f` as bounded composition tracers. At the start of `thermal_step`, the code copies `Hvol` into `Hadv` and inverts that enthalpy using the currently advected material fractions. In the evaporation-enabled path, the P1P source then reconstructs an interface temperature and queries the Ag property table before entering the thermal nonlinear solve.

This arrangement creates a concrete audit question: do enthalpy and composition remain locally compatible after their respective advection updates, especially near discontinuous material boundaries and nearly empty VOF cells? It is not yet evidence that the fluxes are inconsistent. The no-evap smoke failure occurs inside the later nonlinear heat solve; the C stop occurs earlier during P1P property validation. Logs currently do not identify a shared failing cell or show that the same numerical defect triggers both events.

For C, the replayed rejected cells have tiny stored material fractions (`f≤1.17×10⁻⁹`, `Ag≤5.06×10⁻¹¹`, `Ti₃SiC₂≤1.12×10⁻⁹`). With air volumetric heat capacity `1206 J m⁻³ K⁻¹`, their logged `Hadv` maps to a gas-dominated estimate of roughly `4293–4911 K` before the interface-gradient extrapolation. This explains why those cells exceed the property table once classified as exposed Ag interfaces, but it does not establish how their enthalpy arrived there.

For the no-evap smoke, the archived first-step history shows close global energy closure and mass conservation. The second-step log records only the nonlinear update magnitude and stop time; it does not contain per-cell residuals, cellwise `Hadv`, composition, or the iteration-by-iteration Newton history needed to diagnose the convergence failure. Therefore the root cause remains unresolved.

### Offline read of the preserved no-evap failed dump

The archived 99 MiB failed dump was decompressed and read without calling `run()` or advancing the solution. All `311,296` cells restored. Using the archived 300–2300 K Ti₃SiC₂ property table, the local `Hadv`/composition inversion spans `1239.29–1409.79 K`; there are no cells above 2000 K, and the stored component identity `Ag + Ti₃SiC₂ + pore = f` has zero maximum error at stored precision. The failed dump therefore does not show the 5000 K sentinel or the near-empty gas/Ag-trace condition observed in C.

The saved terminal nonlinear iterate has `Tfield=1227.82–1409.83 K` and `Tguess=1178.51–1409.83 K`. Its largest `|T(Hadv,composition)−Tfield|` is `59.793 K`, at a full Ti₃SiC₂ cell near `(x,y,z)=(18.8661,10.4141,0.1509) nm`. The maximum 24-bisection enthalpy inversion residual is `15,820.32 J m⁻³` (`4.31×10⁻⁶` of the local `Hadv`). These results localize the failed endpoint to nonlinear thermal iteration rather than an over-range mixture inversion in the stored state, but they do not identify why the Newton/Poisson iteration failed.

At the terminal failed iteration, `Hvol` is the enthalpy of the current nonlinear candidate while `Hadv` retains the advected target enthalpy. Their maximum difference (`1.2896×10⁹ J m⁻³`) is therefore not an advection-conservation residual and must not be used as proof of enthalpy/composition transport mismatch. The next diagnostic needs iteration-by-iteration residual and cellwise field logging on a separate diagnostic source copy. Summary values are in `P1P_NO_EVAP_FAILED_DUMP_FIELD_AUDIT.csv`; the offline reader is `inspect_no_evap_failed_dump.c` with `run_no_evap_failed_dump_audit.sh`.

## Required next diagnostic, still without production changes

1. Preserve the no-evap smoke logs and failed dump as the regression fixture; do not rerun a heating campaign.
2. Inspect the archived no-evap failed dump with an offline diagnostic build. At the failure state, report the extrema and location of `Hvol/Hadv`, `ag_marker`, `sc_marker`, `pore_marker`, `f`, reconstructed temperature, and local material simplex error.
3. Add diagnostic-only logging in a separate source copy to distinguish: after tracer advection; after `Hadv`/composition temperature inversion; and each thermal nonlinear iteration. Record cellwise and global nonlinear residuals, not only the max temperature update.
4. Audit whether the two transport mechanisms preserve the same composite-volume and enthalpy flux in a minimal advection fixture. Then test enthalpy inversion round trips for gas-only, Ag-only, skeleton-only, and mixed cells.
5. Re-evaluate the archived C cells against that baseline audit. Keep the C property guard as a separate failure until evidence connects it to the baseline nonlinear failure.
6. Do not change thresholds, clip the enthalpy, extrapolate properties, modify the production solver, or resume formal flow until the baseline behavior is understood and the planned numerical checks pass.

## Reproducibility references

- Handoff requirement to resolve the evaporation-off second-step failure first: `../../source_data/user_implementation_pack/P1P_Ag_evaporation_implementation_pack_IMPLEMENTATION_HANDOFF.md`, section 9.
- Source audit statement: `../../source_data/user_implementation_pack/P1P_Ag_evaporation_implementation_pack_SOURCE_AUDIT.md`, “Main implementation blockers”.
- Smoke evidence and run limitations: `../../tests/solver_smoke/README.md` and preserved stderr/history files under `../../tests/solver_smoke/`.
- C-case cell replay: `../p1p_cell_failure_diag_20261107/P1P_C_FAILURE_REPLAY.md`, `debug.stderr`, and `P1P_C_FAILURE_CELL_DIAGNOSTIC.csv`.
- Offline no-evap failed-dump field summary: `P1P_NO_EVAP_FAILED_DUMP_FIELD_AUDIT.csv`; source and commands: `inspect_no_evap_failed_dump.c` and `run_no_evap_failed_dump_audit.sh`.
- Production source SHA256 at reviewed commit: `098a64e3d7123114929d1cd09f27ea7ea1d96d2eac6b9c4091f7b7afbcd22bdd`.
- Archived no-evap failed-dump SHA256 (compressed): `72315ca41dc8c54abadf3a5869768d8f1daa460516443dc1650d78dcd607c5a0`; diagnostic source SHA256: `a6e199c29a503dfdc1118fc82d163477f4e4018ff2b61d59bfdafaf6d809237f`; property table SHA256: `0771fc059d3a597fdf524e3e30cdb54bb47a8850e7e8ffb158852e523645e67c`.
