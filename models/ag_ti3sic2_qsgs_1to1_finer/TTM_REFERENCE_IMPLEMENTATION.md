# Ag-only TTM reference implementation

Status: numerical reference and unit-tested operator prototype. It is **not yet a LAMMPS fix or a production MD input**. No composite MD or heating trajectory has been run with this code.

## What it implements

- The grid covers only the material thickness `z=0–158.6115 Å`; vacuum is excluded. x and y are periodic. Electron flux across both z faces is zero, which represents the open free surface without wrapping electrons through vacuum.
- Only Newton-layer Ag (`type 4`, `z≥40 Å`) contributes to the electron fraction, source, and electron–lattice transfer. Ti₃SiC₂ receives no direct source or Ag electronic parameters.
- Atom volume and source energy are transferred to the mesh with cloud-in-cell weights. Occupancy overfills are locally redistributed while preserving total Ag volume and source power.
- The electronic state is stored as tabulated `Ue(Te)`. Ce, G, ke, and Ue are read from the byte-verified candidate CSV. The table is interpolated without high-temperature extrapolation; values above 20,000 K stop with an error.
- Electron diffusion uses finite-volume face fluxes with a harmonic mean of `fAg·ke`. Equal and opposite face fluxes conserve electron energy; empty Ag cells have no electronic heat capacity or conduction.
- The source is the confirmed lateral Gaussian × depth exponential, multiplied by the Ag-associated cell fraction. It is injected into the electron energy field, not directly into atom velocities.
- The mean-field electron–lattice term is `G(Te)·fAg·Vcell·(Te−Tl)`. The returned per-cell energy is applied to Newton-layer Ag kinetic energy through deterministic velocity scaling. The implementation checks the global electron-plus-lattice energy ledger.
- A separate controller turns heat off after a connected liquid cluster persists for 3 ps and the fitted pool-volume drift stays within 5% over 2 ps. It immediately turns heat off and stops the run if the skeleton disconnects, damage exceeds 5%, surface recession exceeds 3 Å, or any atom escapes. These are provisional operating thresholds, not material constants.

## Low-temperature assumption

The supplied Ag table begins at 300 K. To permit local kinetic-temperature fluctuations below 300 K without clamping electron energy, the reference solver has an **opt-in** continuation from 1–300 K: `Ce(T)=Ce(300)·T/300`, `G(T)=G(300)`, `ke(T)=ke(300)`, with `Ue` integrated consistently. Only the linear Ce form follows the Sommerfeld low-temperature form; constant G and ke are explicit convenience assumptions. This continuation is not measured by the uploaded data and must be included in any parameter-sensitivity study. The original CSV remains unchanged.

## Limitations before any heating run

1. The deterministic mean-field transfer conserves energy but omits Langevin fluctuations. It is a TTM-MD operator candidate, not a drop-in equivalent to `fix ttm`.
2. The 300 K table values are cross-paper fitted candidates, not a self-consistent validated parameter set.
3. The controller consumes liquid-cluster size, pool volume, skeleton connectivity/damage, surface recession, and escaped-atom metrics. The CNA-plus-1-ps-displacement liquid flag has been checked against repository-EAM pure Ag solid/liquid references (see [classifier calibration](pure_ag_classifier/README.md)); the classifier and cluster metric still need validation on composite trajectories, and the online trajectory-to-controller monitor is not yet implemented.
4. A one-step 10 Å geometry smoke test at 0.5 fs required 348 explicit electron substeps. This proves only numerical operation and energy accounting; it also shows that a Python reference solver is too costly for a 100 ps production trajectory. Production needs a compiled LAMMPS fix or a validated implicit/accelerated implementation.
5. The composite cold baseline and current calibration-v0 interface potential have not passed force/equilibrium checks. The Ti₃SiC₂ potential is also not validated for weld-like thermal cycling. The intended observable is a microscopic local Ag melt pool, confinement of molten Ag by the continuous Ti₃SiC₂ scaffold, and local resolidification after heat-off, not bulk melting. The composite trajectory can only be interpreted conditionally on these potentials; pure-Ag calibration does not validate the interface binding or backbone response.

## Reproduction

From this directory:

```bash
python -m unittest -v test_ag_ttm_reference.py test_melt_pool_control.py
python preflight_ttm_grid.py
python ttm_geometry_one_step_smoke.py
```

The last command uses the real atom coordinates and an explicitly seeded synthetic 300 K velocity field. It runs no LAMMPS command and takes no MD step.
