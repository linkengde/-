# Ag phase coverage

## Confirmed DFT-labelled coverage

The archive contains seven unique 32-atom FCC Ag control geometries: one baseline, ±2% volume, ±0.5% shear and two displacement variants. Their nearest Ag–Ag distance range is 2.873–2.911 Å. Multiple k meshes for one fixed geometry are deduplicated. These controls support basic solid-state and local displacement method tests; they do not by themselves define a temperature-dependent Ag model.

## Missing coverage

No converged DFT-labelled liquid Ag configurations were found. There is no temperature/density-spanning liquid set, liquid-solid coexistence interface, or DFT-labelled freezing/solidification sequence. The archived 77-atom gap structures are geometry-only and are not Ag phase labels. MACE, SevenNet or MatterSim predictions cannot be used as DFT labels.

## Decision

`AG_SOLID_DATA_READY = INSUFFICIENT_DATA`: a small solid control set exists but its independent-parent and thermodynamic-response coverage is limited.

`AG_LIQUID_DATA_READY = FAIL`: liquid density/temperature and liquid-solid environments have zero confirmed DFT-labelled geometry families. Additional DFT is required for a model intended to melt and resolidify Ag.
