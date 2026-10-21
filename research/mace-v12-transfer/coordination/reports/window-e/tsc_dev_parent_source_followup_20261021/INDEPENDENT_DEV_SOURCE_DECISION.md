# Independent DEV source follow-up

Date: 2026-10-21
Owner: E
Scope: non-TEST source and geometry-family qualification; no new DFT, relaxation, fitting, role reassignment, or TEST access in this follow-up.

## Decision

`INDEPENDENT_CRYSTALLOGRAPHIC_DEV_PARENT=NOT_AVAILABLE_FROM_AUDITED_SOURCES`.

The C-window crystallography audit compared the 1967 Jeitschko–Nowotny, 1998 COD/Kisi, and 2000 Rawn structures. They share the same Ti3SiC2 `P6_3/mmc` prototype, Wyckoff sites, composition, and layer sequence. The Rawn geometry is nearly coincident with COD. The Jeitschko reconstruction differs more in the C internal coordinate, but it is still the same crystallographic parent family. Neither reconstructed CIF contains DFT energy/force labels. Running new DFT on these CIFs would add same-family geometry data, not an independent parent-family DEV label.

The 2026 Ag/Ti3SiC2 paper supplies scalar interface-spacing/adhesion values, but not complete coordinates or per-atom forces; it cannot supply a reproducible DFT force-label parent.

## Existing E prismatic pair

The E prismatic `(10-10)` slab is a different surface-orientation geometry family from the basal slabs, but it was built from COD 9009647. Its 96-atom fixed-geometry reference and terminal-Ti displacement labels have already been completed and archived. The same-geometry TSCBOP `run 0` comparison is also complete. These are useful **geometry-family diagnostics** and can support a facet-held-out comparison if a future model evaluation protocol explicitly freezes that role. They do not meet a strict independent crystallographic-source blind-test claim. Existing TRAIN/DEV/TEST roles remain unchanged.

The archived pair has:

- REF input SHA-256: `60096a122270e57ac2fb1cbbc2c0d1b087eab12127bcee087f30d244dc8900e1`; geometry SHA-256: `30a4b161cd9a912c90006f3c0d35978e4d417306734f3fc982c9c8684fa600d3`.
- Ti displacement input SHA-256: `710b0abdd20f0d59103ebfe061ad65c834bae7a4b2f5a300453ba1527e63e030`; geometry SHA-256: `33aa1032f0ec404b5ca46976f82da5c4e1dc29cff8cebb76cf3f25d8a7a92e49`.
- Both structures have 96 atoms (`Ti48 Si16 C32`) and share source family `surface_prismatic_10m10_cod9009647_8layer`.
- The reference TSCBOP absolute force-vector RMSE is 1.9152 eV/A; the displaced point is 1.9130 eV/A. The paired differential-response RMSE is 0.01348 eV/A over the full slab and 0.04101 eV/A within 3.3 A of the displaced atoms. This one response direction does not establish broad transferability and does not remove the large absolute Ti/C force bias.

## Recommended next step

Do not spend another DFT job on the Rawn or Jeitschko reconstruction solely to claim independent validation. Keep the completed prismatic pair as a geometry-family diagnostic candidate. If strict independent-parent validation remains a gate, obtain a genuinely new parent with complete structure provenance and freeze it before any model scoring; until then mark independent validation unavailable. A grouped holdout by physical environment (entire facet/termination or independently generated configuration family) is a useful secondary diagnostic, but it must be labeled as such and must not be described as an independent crystallographic-source blind test.

No DFT, MD, TTM, parameter optimization, or potential-file modification was performed for this follow-up. No data role was reassigned. This follow-up did not access TEST inputs or labels; the earlier access disclosure remains in the prior E diagnostic report and is not reinterpreted here.

## Evidence sources

- C audit commit: `e999abccf69f717f29d255377b006d4d4a258266`, report `v19_three_paper_crystallography_interface_audit_20261011/REPORT.md` and its `SHA256SUMS.txt`.
- E geometry-family DFT/TSCBOP audit commit: `7df8ac8daa3811e755985d69e5a5fceffe4889e9`.
- E TSCBOP static comparison commit: `c02631444b5d41a71b58307aed83b0fb421929a4`.
