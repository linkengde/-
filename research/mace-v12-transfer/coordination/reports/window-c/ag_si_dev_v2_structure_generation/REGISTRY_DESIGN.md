# Ag/Ti₃SiC₂(0001)-Si registry candidates

## Literature basis and matching decision

The design-method reference is A. H. M. Abdul Wasey et al., “First principles electronic structure of coincidence site epitaxial Ag/Si(111) interface,” *physica status solidi (b)* **250**(7), 1313–1319 (2013), DOI `10.1002/pssb.201248542`. The user-supplied primary PDF reports the fcc Ag cubic lattice parameter `a₀=4.09 Å`; thus ideal Ag(111) nearest-neighbor spacing is `a₀/√2=2.892067 Å`.

That paper's **4×Ag(111)/3×Si(111)** match is specific to diamond Si(111): using its `a(Si)=5.43 Å`, the Si(111) surface primitive length is `5.43/√2=3.8396 Å`; `3×3.8396=11.5188 Å`, close to `4×2.8921=11.5683 Å` (about 0.43% mismatch using rounded reported inputs). It is not transferable as the same 4:3 ratio to Ti₃SiC₂(0001): for `a=3.071 Å`, `3×a=9.213 Å`, about 25.6% smaller than `4×2.8921 Å`. Applying that exact ratio here would impose a grossly mismatched cell.

For the target MAX-phase interface, this package therefore uses the target-specific **1×1 coherent hexagonal cell** described in the reviewed direct interface study by Zhang et al., *Coatings* **16** (2026) 995, DOI `10.3390/coatings16080995`: `a=b=3.071 Å`, `γ=120°`, Ag coherently matched to the substrate, seven Ag(111) atomic layers. The COD/Kisi fractional substrate coordinates are reconstructed in this cell. The 3.071 Å value and seven-layer construction are literature-derived; the atom-coordinate reconstruction and registry variants below are generated here, not copied from an author structure file.

Coherently matching the unstrained Ag nearest-neighbor spacing 2.892067 Å to 3.071 Å gives +6.187% in-plane strain. To keep the unrelaxed Ag ABC slab internally ideal fcc, the generator scales the ideal Ag cubic parameter isotropically to `√2×3.071=4.343050 Å`, giving Ag(111) interlayer spacing `2.507461 Å`. This is an explicit starting-geometry choice, not an optimized strain response.

## Registry definitions

Translations are fractional coordinates in the common γ=120° substrate vectors. The terminating Si site is at `(0,0)`. One Ag atom occupies each layer of the 1×1 cell. The seven-layer Ag stacking is ABCABC A, with incremental in-plane offsets `(0,0)`, `(2/3,1/3)`, `(1/3,2/3)` repeated; the listed registry shift is added to every Ag layer.

| ID | Site | Bottom Ag shift `(u,v)` | Geometric definition |
|---|---|---:|---|
| R1 | top / OT-like | `(0,0)` | Bottom Ag directly above terminating Si |
| R2 | bridge / MT-like | `(1/2,0)` | Bottom Ag above the midpoint of neighboring surface Si sites along `a₁` |
| R3 | threefold hollow / HCP-like | `(2/3,1/3)` | Bottom Ag above the immediate subsurface Ti lateral site |
| R4 | alternate hollow / FCC-like | `(1/3,2/3)` | Bottom Ag above the next subsurface C lateral site |

The direct interface study names top (OT), bridge (MT), and HCP registry classes. R4 is an additional deterministic alternate hollow constructed for this scan; it is not claimed as a fourth registry reported by either paper. “HCP-like” and “FCC-like” describe the registry relative to the reconstructed MAX-phase subsurface sequence, not an assertion that the substrate itself is elemental hcp/fcc.

## Vertical construction and screening

The mean-plane gap is fixed at `2.600 Å` for all four unrelaxed geometries, measured between the mean z of the outermost Si termination plane and the mean z of the bottom Ag plane. This single common value makes registry geometry comparisons well-defined. It is **not** a fitted equilibrium spacing or a conversion of the 2013 paper's local 2.53–2.72 Å Ag–Si bond lengths. The direct interface paper's tabulated `d₀` values also are not reused as mean-plane gaps because its coordinate-level convention cannot be mapped from the available structure data.

Minimum pair distances are computed from the written extended XYZ using ASE 3.29.0 minimum-image distances in x/y only (`PBC=T,T,F`). The acceptance checks are `min(Ag–Si) ≥ 2.0 Å` and no Ag–Ti or Ag–C pair below the same severe-overlap floor. The latter is only a gross-overlap screen; passing it does not establish favorable bonding, relaxed stability, or physical accuracy. All four candidates pass this geometric check; exact minima are in `REGISTRY_SCREENING.csv`.

## Data role and scope

The substrate parent is COD 9009647/Kisi 1998, which shares the same crystallographic parent family as previously scored COD-derived Ag–Si structures. Although these are new 20-atom coordinates and include a different 1×1/seven-layer construction, the family is not source-independent. The structures remain `GEOMETRY_ONLY_NOT_ASSIGNED_TRAIN_DEV_TEST`; this report does not promote them to independent DEV data or alter any existing role. A future DFT calculation requires separate authorization. No DFT, optimization, MD, MACE inference, or training was performed here.

**Geometry status:** `READY_FOR_DFT_SCREENING`.

**Independent DEV status:** `NOT_ESTABLISHED`.
