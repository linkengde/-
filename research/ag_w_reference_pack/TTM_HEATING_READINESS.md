# Ag-W TTM heating readiness

## Requested geometry and heating method

The candidate model uses x/y periodic boundaries and a free z surface (`p p f`), with fixed atoms at 0–20 Å, a thermostat layer at 20–40 Å, a Newton layer at 40–200 Å, and empty sputter/evaporation monitor regions at 200–250 Å and 250–300 Å. Laser heating is to be done with TTM; no direct thermostat pulse is substituted.

## Why there is no TTM trajectory yet

The installed LAMMPS build contains TTM styles, but the existing project preflight for `fix ttm/mod` found that its boundary assumptions reject nonperiodic z. Changing z to periodic would make atoms that leave the free surface wrap through the vacuum to the bottom, so that stock fix cannot be used with this geometry as-is.

There is an Ag-only finite-volume TTM reference solver in `models/ag_ti3sic2_qsgs_1to1_finer/ag_ttm_reference.py`. It supports periodic x/y and zero-flux open z faces and has been coupled to short LAMMPS Python API pilots. It is a custom numerical prototype, not a compiled LAMMPS fix, and it does not provide W parameters or a validated Ag–W mixture closure. It is a possible implementation starting point; the Ag-only table from that separate project must not be presented as an Ag–W parameter set.

The Ag–W reference pack currently has no verified temperature-dependent electronic heat capacity `Ce(Te)`, electronic thermal conductivity `ke(Te)`, or electron–lattice coupling `G(Te)` tables for W, nor an approved rule for electron transport and coupling in cells containing both phases. The laser source also has no Ag–W-specific fluence or absorbed energy, pulse duration, lateral profile, or optical penetration depth specified. These inputs are needed before the TTM prototype can be adapted and its energy ledger checked.

For those reasons, the new QSGS geometry report records TTM as requested but `NOT_RUN`. A zero-source 300 K exploration is also pending until the raw close contacts are relaxed; the initial perfect-lattice assignment has a 0.08 Å Ag–W distance and must not be integrated directly.

## TTM acceptance inputs

Before a small laser-heating pilot, record:

1. Literature-backed `Ce(Te)`, `ke(Te)`, and `G(Te)` for Ag and W, with units and temperature ranges.
2. A coupling implementation that supports `p p f`, the selected atomistic potential, the Ag/W type distinction, and a free surface. The existing Ag-only finite-volume solver can inform this implementation but does not meet the Ag–W requirements yet.
3. Absorbed fluence or total deposited energy, pulse duration, lateral source profile, absorption depth, and whether the 20–40 Å layer is a heat sink during the pulse.
4. Energy-accounting checks for source deposition and equal-and-opposite electron-lattice transfer before any trajectory is interpreted.

The current structure is only an unrelaxed geometry candidate. It has not been energy minimized, equilibrated at 300 K, heated, or used for production MD.
