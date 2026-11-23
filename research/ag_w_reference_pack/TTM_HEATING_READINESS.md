# Ag-W TTM heating readiness

## Requested geometry and heating method

The candidate model uses x/y periodic boundaries and a free z surface (`p p f`), with fixed atoms at 0–20 Å, a thermostat layer at 20–40 Å, a Newton layer at 40–200 Å, and empty sputter/evaporation monitor regions at 200–250 Å and 250–300 Å. Laser heating is to be done with TTM; no direct thermostat pulse is substituted.

## Why there is no TTM trajectory yet

The installed LAMMPS build contains TTM styles, but the existing project preflight for `fix ttm/mod` found that its boundary assumptions reject nonperiodic z. Changing z to periodic would make atoms that leave the free surface wrap through the vacuum to the bottom, so that stock fix cannot be used with this geometry as-is.

There is an Ag-only finite-volume TTM reference solver in `models/ag_ti3sic2_qsgs_1to1_finer/ag_ttm_reference.py`. It supports periodic x/y and zero-flux open z faces and has been coupled to short LAMMPS Python API pilots. It is a custom numerical prototype, not a compiled LAMMPS fix, and it does not provide W parameters or a validated Ag–W mixture closure. It is a possible implementation starting point; the Ag-only table from that separate project must not be presented as an Ag–W parameter set.

The spatial heat-source formula and candidate values were migrated into an independent Ag–W static preflight. Its 10 Å Ag mask was built from the original, overlapping type-2 atom coordinates; at `S(t)=1`, the mapped Newton-layer Ag power was 9,784.984 eV/ps. This is a geometry-specific on-state power, not a pulse energy, and is **stale for any future reconstructed coordinates**. Do not reuse this power or mask: recalculate after a candidate passes static structural acceptance. The historical report, JSON, and mask/source arrays are in `models/ag_w_qsgs_200A/ttm_source_preflight/`. The legacy Ag–Ti₃SiC₂ files were not changed.

The Ag–W reference pack now contains public XTANT-3 elemental Ag and W tables for `Ce(Te)`, `G(Te)`, and `ke(Te)` in `ttm_data/XTANT3_TTM_parameters_v3/`. However, the W `Ce` and `G` tables span 500–19000 K, the W `ke` table spans 300–50000 K, and the dataset warns that room-temperature values may not always be reliable. The imported `G` and `ke` tables are indexed by `Te` only, not full `(Te,Tl)` surfaces. A W `ke(Te,Tl)` Drude-style closure is available in Daraszewicz et al. (2014) as a model, and the W low-temperature `G_0` values remain method-dependent (see `literature_20261123/README.md`). There is still no approved mixed-cell Ag/W transport rule or Ag/W electronic interface conductance. The exploratory spatial source profile is specified, but `S(t)` has no pulse width or temporal envelope, and no Ag–W-specific absorbed fluence is specified. These inputs are needed before the TTM prototype can be adapted and its energy ledger checked.

For those reasons, TTM is `NOT_RUN`. V3 local slabs are mechanically relaxed, but v1/v2 did not produce a valid full-size starting structure; the raw perfect-lattice assignment has a 0.08 Å Ag–W distance and must not be integrated directly. Qualitative use of the potential is accepted for graduate-level exploration, but any 300 K MD still requires a full 3D candidate to pass static structure, density, and connectivity checks. TTM also requires the missing W data and mixed-phase closure. See `models/ag_w_qsgs_200A/interface_reconstruction_v3/LOCAL_INTERFACE_V3_AUDIT.md`.

## TTM acceptance inputs

Before a small laser-heating pilot, record:

1. Select and reconcile the imported literature-backed `Ce(Te)`, `ke(Te)`, and `G(Te)` branches for Ag and W, with units and temperature ranges; preserve the model disagreement as sensitivity branches.
2. A coupling implementation that supports `p p f`, the selected atomistic potential, the Ag/W type distinction, and a free surface. The existing Ag-only finite-volume solver can inform this implementation but does not meet the Ag–W requirements yet.
3. Temporal envelope `S(t)` or absorbed fluence/total deposited energy, pulse duration, confirmation that the candidate lateral radius/depth are appropriate for Ag–W, and whether the 20–40 Å layer is a heat sink during the pulse.
4. An Ag/W interface electronic conductance sensitivity model, since no direct Ag/W value was identified; energy-accounting checks for source deposition and equal-and-opposite electron-lattice transfer before any trajectory is interpreted.

No full-size 3D structure has passed static acceptance. The original and v1/v2 candidates have not been equilibrated at 300 K, heated, or used for production MD. V3's small local slabs were zero-temperature EFS minimizations only.
