# Ag–W v6 TTM link preflight

## Decision

The synthetic Ag–W TTM grid smoke check passes. The finite-size atom-to-grid preflight fails, so the LAMMPS/TTM source and matched no-source trajectories were **not** run. The phase-volume check was not relaxed or clipped.

This is a software/mapping blocker, not evidence that the TTM energy solver failed. No laser heating, melting, or cooling result was generated.

## Cold state used

- v6 fixed-ID 1.0 ps NVT plus 0.05 ps NVE continuation.
- 482,601 atoms retained; original 46,117 bottom IDs remained fixed.
- NVT endpoint mobile temperature: 293.79 K (Ag 289.35 K; W 301.08 K).
- NVE temperature change: −0.084 K; total-energy drift: +0.00038 eV.
- NVT temperature was still rising by 8.51 K over the final 0.1 ps; this remains an exploratory candidate, not a fully equilibrated or static-pass state.
- Exact start data hash and restart-dump hash are recorded in `interface_reconstruction_v6_exploratory/exploratory_ttm_state_manifest.json`.

## Checks completed

1. `run_agw_ttm_grid_smoke.py` passed on its synthetic 8×8×4 two-phase grid. Ag-only source mapping matched its reference; unlike-phase interface energy transfer was equal and opposite; combined electron/lattice energy residual was `−1.48e-14 eV`.
2. The LAMMPS pilot loaded the v5 EFS data and the v6 NVE endpoint restart, preserving the same fixed IDs. It completed the 20-step 300 K setup segment.
3. The 10 Å CIC map then rejected the actual geometry because its peak `Ag+W` phase fill was `1.31585`, above the cell capacity of 1.0. No TTM grid was initialized for this state; the requested source/no-source comparison did not begin.
4. A mapping-only mesh scan gave these maximum fills:

| Spacing (Å) | Mesh | Peak phase fill | Result |
| --- | --- | ---: | --- |
| 5×5×5 | 40×40×32 | 1.74680 | Reject |
| 10×10×10 | 20×20×16 | 1.31600 | Reject |
| 20×20×20 | 10×10×8 | 1.08667 | Reject |
| 40×40×40 | 5×5×4 | 1.02725 | Reject |
| 20×20×40 | 10×10×4 | 1.05161 | Reject |
| 40×40×20 | 5×5×8 | 1.04867 | Reject |
| 50×50×40 | 4×4×4 | 1.02421 | Reject |
| 100×100×40 | 2×2×4 | 1.01451 | Reject |
| 100×100×20 | 2×2×8 | 1.02635 | Reject |
| 100×100×80 | 2×2×2 | 1.00666 | Reject |
| 200×200×40 | 1×1×4 | 1.01153 | Reject |
| 200×200×80 | 1×1×2 | 1.00571 | Reject |
| 200×200×160 | 1×1×1 | 0.99600 | Capacity passes, spatially unusable |

The only passing mesh collapses the whole domain into one cell. It cannot resolve the Ag/W interfaces, the 30 Å spot, or spatial transport, so it is not an acceptable workaround.

## Attempt log

The initial attempts exposed and fixed pilot implementation/runtime issues before reaching the TTM calculation:

- The bundled LAMMPS Python library needs its MPICH runtime directory on `LD_LIBRARY_PATH`; `mpi4py` is absent but not needed for the single-process run.
- A potential-path lookup incorrectly assumed the repository root was the Ag–W pack root. It now searches for the potential file itself.
- Printing the Langevin tally from LAMMPS thermo during a setup `run 0` crashed the bundled runtime; the coupled loop now turns off thermo printing and reads the tally directly.
- A NumPy view into LAMMPS atom arrays became stale after NVT reallocation; the pilot now refreshes views after dynamics steps.
- With those corrections, the pilot reached the intended geometry check and stopped cleanly at the overfill guard.

Retained LAMMPS logs:

- `agw_ttm_attempt_thermo_segfault.lammps.log`
- `agw_ttm_attempt_array_view_segfault.lammps.log`
- `agw_ttm_attempt_mapping_overfill.lammps.log`

## Next action

Keep the TTM bridge’s phase-volume guard. Resolve the local atom density/volume mapping for this candidate, then regenerate the Ag mask and rerun the mapping check. Only after the 10 Å (or otherwise justified, spatially resolved) map satisfies the cell-volume constraint should matched no-source and low-integrated-source coupling runs be attempted. The static EFS audit remains `STATIC_EFS_FAIL`; the exploratory manifest does not change that result.
