# Minimum additional DFT plan (planning only)

No calculations are authorized or performed by this audit. First recover metadata and standardize the confirmed existing labels; do not recompute structures just to duplicate a different k mesh.

## Priority sequence

| Priority | Planned independent DFT configurations | Main constraint | Approx. count |
|---|---|---|---:|
| P0 backbone response | Ti₃SiC₂ volume −2%/+2%; shear −0.5%/+0.5%; Si displacement ±0.02 Å from the existing bulk parent | elastic response, Si force response, parameter separation; existing Ti/C displacements do not replace these | 6 |
| P1 clean surface | one independent Si-terminated slab baseline plus 3 small symmetry-distinct under-coordination/thermal perturbations | dangling/low-coordination response and a real surface family; 114/135 interface-derived subsets are not clean surfaces | 4 |
| P2 mixed interface | two distinct registries/terminations, each at compressed, near-equilibrium and separated gap | four-element mixed interface energy/force and range; reuse existing verified pair-family labels where geometry/method match, but do not treat them as liquid contact | up to 6 new (up to 9 if no old point is compatible) |
| P3 Ag liquid | 4 decorrelated snapshots spanning at least two target liquid densities/temperatures | Ag liquid energies and forces; no DFT-labelled liquid snapshot exists | 4 |
| P4 liquid-solid contact/coexistence | 6 decorrelated configurations across two contact gaps and solid/liquid local states, assigned by whole parent trajectory | Ag liquid touching solid backbone and incipient freezing; no such labelled family exists | 6 |

Indicative total: 26 new configurations if 3 of the P2 points can be reused after strict method/geometry audit; otherwise up to 29. A later independent TEST family must be generated and kept sealed; the table's counts include independent families for that purpose only if C assigns a new parent family before computation. Do not split adjacent frames from one trajectory between splits.

## Cost estimate

The B/V19 archive contains 93 unique `summary.json` runtime records with atom count and elapsed time. For 1–48 atoms, median elapsed wall time is about 4,206 s (70 min) and observed maximum about 15,229 s (4.23 h). For 49–80 atoms, median is about 5,334 s (89 min) and observed maximum about 15,766 s (4.38 h). These are historical GPAW/PBE runs with differing k meshes and environments; they are planning references, not guarantees for a new node.

If all 26–29 new points happened to remain within the archived ≤80-atom regime, a serial median-based envelope is roughly 38–43 wall hours and an observed-maximum envelope roughly 114–127 wall hours. Slab/contact structures may exceed 80 atoms, and no validated runtime extrapolation is available above that range. The 135-atom 2×2×1 attempt previously hit a 28 GiB memory guard before SCF, so larger models require a small approved resource pilot and likely smaller cells. These estimates exclude reruns caused by convergence failures and are not an execution authorization.

## Role isolation

C must freeze parent-family roles before any new DFT. TRAIN can contain selected bulk response and one surface/interface family; DEV must be a separate termination/registry and separate liquid-contact parent; TEST must be a newly constructed independent interface or coexistence parent not used in fitting or model selection. Keep all siblings, k-mesh variants and frames from one source family together. Do not reassign existing frozen test data.
