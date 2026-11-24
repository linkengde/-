# Ag-W TTM heating readiness

## Requested geometry and heating method

The candidate model uses x/y periodic boundaries and a free z surface (`p p f`), with fixed atoms at 0–20 Å, a thermostat layer at 20–40 Å, a Newton layer at 40–200 Å, and empty sputter/evaporation monitor regions at 200–250 Å and 250–300 Å. Laser heating is to be done with TTM; no direct thermostat pulse is substituted.

## Why there is no TTM trajectory yet

The installed LAMMPS build contains TTM styles, but the existing project preflight for `fix ttm/mod` found that its boundary assumptions reject nonperiodic z. Changing z to periodic would make atoms that leave the free surface wrap through the vacuum to the bottom, so that stock fix cannot be used with this geometry as-is.

There is an Ag-only finite-volume TTM reference solver in `models/ag_ti3sic2_qsgs_1to1_finer/ag_ttm_reference.py`. It supports periodic x/y and zero-flux open z faces and has been coupled to short LAMMPS Python API pilots. It is a custom numerical prototype, not a compiled LAMMPS fix, and it does not provide W parameters or a validated Ag–W mixture closure. It is a possible implementation starting point; the Ag-only table from that separate project must not be presented as an Ag–W parameter set.

The spatial heat-source formula and candidate values were migrated into an independent Ag–W static preflight. Its 10 Å Ag mask was built from the original, overlapping type-2 atom coordinates; at `S(t)=1`, the mapped Newton-layer Ag power was 9,784.984 eV/ps. This is a geometry-specific on-state power, not a pulse energy, and is **stale for any future reconstructed coordinates**. Do not reuse this power or mask: recalculate after a candidate passes static structural acceptance. The historical report, JSON, and mask/source arrays are in `models/ag_w_qsgs_200A/ttm_source_preflight/`. The legacy Ag–Ti₃SiC₂ files were not changed.

The Ag–W pack now has public XTANT-3 Ag and W `Ce(Te)`, `G(Te)`, and `ke(Te)` tables in `ttm_data/XTANT3_TTM_parameters_v3/`; original values and units are retained. The user-uploaded README and license notices are also archived beside the tables. The standalone `Ce` files are from dynamical simulations at the gamma point, while the `Ce` columns inside K tables are from ideal crystals with a 7×7×7 k-point mesh; these definitions are not interchangeable. At sampled common temperatures, XTANT W `G` is 3.8–17.9 times lower than the Lin et al. digitized curve. The source Zenodo record warns that room-temperature values may not always be reliable. See `ttm_data/XTANT3_AG_W_DATA_AUDIT_20261124.md`, `ttm_data/XTANT3_AG_W_SOURCE_AUDIT.json`, and `literature_20261123/README.md` for the table audit and provenance.

`ttm_data/AgW_TTM_baseline_v1/` now records the selected exploratory baseline: Lin 2008 `Ce/G` for both Ag and W, Daraszewicz 2014's W `ke(Te,Tl)` Drude closure, and Ag's XTANT-3 `kappa_tot(Te)` table. The W low-temperature Lin `G` check is outside Daraszewicz's measured 1σ `G₀` interval; the digitized 300 K Lin point is itself extrapolated, and the baseline keeps this discrepancy visible without rescaling the Lin curve. One Lin-W 0D energy-ledger pilot passes its numerical conservation threshold. The earlier W `Ce/G` source audit remains in `ttm_data/W_EPH_TIMESCALE_AUDIT_20261124.md`; it is a literature audit, not a request to run three TTM trajectories.

The selected baseline is explicit, but it is a practical exploratory closure rather than a fully measured set of `(Te,Tl)` surfaces. Ag conductivity is a `Te`-only XTANT table; W conductivity uses the Daraszewicz Drude expression. Daraszewicz's low-temperature W `G0` anchor does not agree with the Lin digitized 300 K extrapolation within its quoted 1σ interval. No Ag/W-specific mixed-cell transport rule or electronic interface conductance is available. The supplied interface-conductance paper covers Au/Pt and Au/Al, not Ag/W. The spatial source profile is recorded, but `S(t)` has no target pulse width or temporal envelope, and no Ag–W-specific absorbed fluence is set.

The inherited spatial Gaussian profile is preserved. Before a 3D pilot, define `S(t)` and absorbed energy, supply an Ag lattice-conductivity choice, choose an explicit Ag/W electronic interface/mixed-cell rule, and implement a solver compatible with `p p f`. The selected baseline and its single-branch W energy ledger are documented in `ttm_data/AgW_TTM_baseline_v1/`; no three-branch 0D trajectory comparison is needed for this baseline decision. The new literature collection includes Ag/W phase and melt references, but the 1988 phase diagram is calculated, not experimentally measured.

Full spatial TTM and atomistic TTM-MD are `NOT_RUN`; only the single W 0D energy-ledger pilot has been run. V3 local slabs are mechanically relaxed, but the v3 full-size candidate still failed the static-start audit after 10,000 EFS minimization steps, with collapsed contacts and atoms above the material surface. The original perfect-lattice assignment still has a 0.08 Å Ag–W distance and must not be integrated. Qualitative use of the potential is accepted for graduate-level exploration, but the 3D candidate must first be rebuilt and pass static geometry, density, and connectivity checks. See `models/ag_w_qsgs_200A/interface_reconstruction_v3/LOCAL_INTERFACE_V3_AUDIT.md` and `models/ag_w_qsgs_200A/interface_reconstruction_v3/qsgs_candidate/FULL_SIZE_STATIC_AUDIT.md` for the structural results.

## TTM acceptance inputs

Before a small laser-heating pilot, record:

1. Use the recorded AgW TTM baseline v1 tables and functions with their stated ranges and provenance. Keep the W low-temperature `G0` mismatch visible; do not silently rescale the Lin curve. No source averaging is used.
2. A coupling implementation that supports `p p f`, the selected atomistic potential, the Ag/W type distinction, and a free surface. The existing Ag-only finite-volume solver can inform this implementation but does not meet the Ag–W requirements yet.
3. Temporal envelope `S(t)` or absorbed fluence/total deposited energy, pulse duration, confirmation that the candidate lateral radius/depth are appropriate for Ag–W, and whether the 20–40 Å layer is a heat sink during the pulse.
4. An Ag/W interface electronic conductance sensitivity model, since no direct Ag/W value was identified; energy-accounting checks for source deposition and equal-and-opposite electron-lattice transfer before any trajectory is interpreted.

No full-size 3D structure has passed static acceptance. The original and v1/v2 candidates have not been equilibrated at 300 K, heated, or used for production MD. V3's small local slabs were zero-temperature EFS minimizations only.
