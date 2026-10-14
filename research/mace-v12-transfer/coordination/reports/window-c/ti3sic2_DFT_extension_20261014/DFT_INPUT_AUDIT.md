# Stage C DFT input audit

## Reviewed method source

The E6 request is `E6_TO_C_DFT_REQUEST.md` from E branch commit `bdba4115569f74578e2fc2e2bee031cebf3f502d`; its manifest SHA is `c081a2dce8faa1a08f9f6a221658dccfa0d4b712cd97b6e9fb511c3fdcf41c61`. Its `DFT_REQUIREMENTS.md` SHA is `e35120d8d09369d0b31db95bc9c37debcfd102bcb9e72183150494cba51a723e`.

E6 and B's verified fixed-cell baseline/response report specify the method to retain:

| Item | Proposed fixed setting | Evidence / status |
|---|---|---|
| Code | GPAW 26.7.0; ASE 3.29.0; gpaw-data 1.2.1 | Installed versions checked in this environment and match B's archived method. |
| Functional / basis | PBE, PAW, plane-wave mode, 500 eV | Match B method manifest SHA `702df71e586a61cd5785afc62653d7135a707032c91b72842b035b62198fd709`. |
| Occupations | Fermi–Dirac, 0.10 eV | Match archived method. |
| SCF thresholds | energy, density and eigenstates: `1e-5` | Match archived method; no threshold is relaxed. |
| Mixer | Pulay β=0.05, `nmaxold=8`, weight=100 | Preserve the reviewed GPAW setup used in the project. |
| SCF cap | `maxiter=160` | E6 method candidate; do not change without a new method audit. |
| Energy labels | Save native extrapolated energy and force-consistent free energy separately | Never substitute one field for the other. |
| Forces / stress | Save complete per-atom forces. Save full stress for bulk points. For slabs, record in-plane components only under a declared slab convention; do not use vacuum-diluted 3D stress as bulk stress. | Required for interpreting the proposed objectives. |
| Geometry | Fixed-geometry single point; no ionic or cell relaxation in this package | No DFT run occurred. |
| Bulk k mesh | `10×10×4` | Match B's baseline and four Ti/C displacement calculations. |
| Slab k mesh | `10×10×1`, GPAW integer Monkhorst–Pack mesh with no added offset | Same lateral cell as bulk, one point along the 20 Å vacuum direction. This is a proposed matched mesh, not newly converged here. A dedicated convergence check is still required before release. |
| D-family k mesh | TBD after an actual independent cell/source is reviewed | Cannot freeze a mesh for missing coordinates. |

## PAW setup identity

Present-element setup files were found in the active GPAW data directory and matched the project hashes:

| Setup | SHA-256 |
|---|---|
| Ti.PBE.gz | `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f` |
| Si.PBE.gz | `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0` |
| C.PBE.gz | `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983` |

Ag setup is not needed for this pure Ti₃SiC₂ dataset.

## Geometry and input completeness

Ten source-backed geometry files were written to `inputs/` and reread with ASE 3.29.0. All have 48 or 100 atoms as listed, Ti/Si/C only, finite coordinates/cells, TTT PBC, unique local IDs, and minimum periodic distance ≥2.074 Å (thus none has a <1.5 Å contact). The four A files are byte-preserved proposal geometries; the B and C files include parent/local-ID mappings and deterministic distortion metadata. Detailed per-file values and SHA-256 are in `DFT_STRUCTURE_LIST.csv` and `geometry_validation.json`.

The already converged parent plus Ti/C ±0.02 Å force-response labels referenced by E6 are existing data and are not counted as new runs in the 14-point list. They are not copied into this package and their current numerical-control roles remain unchanged.

The volume/shear labels add energy, force and stress curvature information. They may improve sensitivity to radial/angle-like parameter blocks, but they do not establish that `lambda2_Ti-C` is independently identifiable. The force/stress response Jacobian and parameter conditioning must be recomputed only after labels are available; E6's present Jacobian is not proof of future identifiability.

The Ti- and Si-terminated slabs are symmetric, unrelaxed cuts with 20 Å vacuum. The terminations have different stoichiometries, so comparing their raw total energies does not measure relative surface stability. A thermodynamic surface-energy comparison would require a separately defined chemical-potential/reference convention and potentially relaxed/reconstructed surfaces; this 4-point proposal only supplies fixed-geometry energy/force response candidates.

No GPAW input script, run manifest or executable DFT authorization is created for the four D rows because the independent structures do not exist in the reviewed source bundle. No SCF was launched for any of the ten geometries, so their force/energy/stress labels are not present. Surface k-point convergence and a resource-approved slab pilot are also outstanding.

## Source and role boundary

The training proposal uses one structural ancestry. The four D points are not generated until separate DEV and evaluator-controlled TEST parent sources are established. Existing numerical-control labels, E archives and all sealed evaluator data remain untouched. See `TRAIN_DEV_TEST_SPLIT.md`.

