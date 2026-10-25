# Ag-only TTM reference implementation

Status: numerical reference and unit-tested operator prototype. It is **not a compiled LAMMPS fix or a production MD input**. A bounded 0.05 ps Model B pilot and a separate 0.05 ps pure-Ag reference-coupling pilot have now been run through the LAMMPS Python API. Neither is a melt-pool run. The new pure-Ag pilot closes its source-corrected combined-energy ledger against matched NVE drift at two mesh resolutions; see [the pilot report](pure_ag_classifier/TTM_MD_PILOT.md). The archived B pilot did not form a stable melt pool and does not validate the primary Morse-interface model; see [the B pilot report](potential_validation/geometry_constraint_B/B_TTM_MD_PILOT_0P05PS.md).

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
4. The one-step 10 Å synthetic-velocity smoke test at 0.5 fs required 348 explicit electron substeps. The 0.05 ps Model B API pilot used actual LAMMPS MD steps and required 242–351 electron substeps per step. After vectorizing lattice-energy handoff and omitting fixed-fixed ceramic calculations under the rigid-obstacle assumption, that pilot took 64.6 s; a simple linear estimate is about 36 h for 100 ps on the same four-thread setup. This is only a planning estimate. A 0.25 ps extension ended with a source-weighted local Ag lattice temperature of 253 K and did not show a melt pool; its corrected energy residual (0.05890 eV) differed from the matched no-source NVE drift (0.05795 eV) by 0.00095 eV. A separate three-times-source 0.25 ps sensitivity reached 2585 K in the electron grid but only 255 K source-weighted lattice temperature. High-temperature performance and long-run accuracy remain unvalidated.
5. The primary Morse-interface composite cold baseline has not passed force/equilibrium checks. Model B is a separate rigid geometric-obstacle exploration using repulsive WCA cross terms. Its bounded 0.25 ps NVE gate failed, and the single permitted 1 ps NVT preparation attempt also failed: maximum Ag force reached 7.45 eV/Å and potential energy continued changing strongly. The planned 1 ps heating pilot has not been run; the current B initial state is rejected for a melting interpretation, with no further relaxation cycle started. The Ti₃SiC₂ potential is not validated for weld-like thermal cycling. The intended observable is a microscopic local Ag melt pool, confinement of molten Ag by the continuous Ti₃SiC₂ scaffold, and local resolidification after heat-off, not bulk melting. Neither the B pilot nor pure-Ag calibration validates real interface binding, ceramic response, or the primary composite's melting behavior.
6. The new pure-Ag sourced reference pilot used 100 actual LAMMPS steps (0.05 ps), the existing Gaussian/exponential source, 0.5 fs MD steps, and 10 Å and 5 Å nominal meshes. After subtracting 260.17268 eV of source input, the combined energy residual was 8.34×10⁻⁵ eV and 8.32×10⁻⁵ eV, respectively; matched no-TTM NVE drift was 5.29×10⁻⁵ eV. This validates only short-time numerical energy handoff for a dense Ag slab under the reference algorithm. It does not validate melting, the QSGS composite, or the candidate material parameters. The 5 Å grid required 33,450 electron substeps for 100 MD steps versus 10,355 at 10 Å, highlighting the reference solver's runtime cost. Full configuration and metrics are in [the pilot report](pure_ag_classifier/TTM_MD_PILOT.md).

## Reproduction

From this directory:

```bash
python -m unittest -v test_ag_ttm_reference.py test_melt_pool_control.py
python preflight_ttm_grid.py
python ttm_geometry_one_step_smoke.py
```

The last command uses the real atom coordinates and an explicitly seeded synthetic 300 K velocity field. It runs no LAMMPS command and takes no MD step.
