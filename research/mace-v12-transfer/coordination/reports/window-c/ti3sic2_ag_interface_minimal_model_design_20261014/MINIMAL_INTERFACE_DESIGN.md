# Minimal Ag/Ti₃SiC₂ Interface Model Design

## Scope and release state

This is a geometry design, not a calculation release. The original 15 fixed DFT inputs remain unchanged and archived in `ti3sic2_ag_interface_fit_dft_inputs_20261014/`.

```text
DFT_STARTED=NO
TSCBOP_MODIFIED=NO
FITTING_STARTED=NO
EXISTING_15_INPUTS_MODIFIED=NO
PROPOSED_COORDINATES_GENERATED=NO
```

The design retains a real Ag(111)/Si-terminated Ti₃SiC₂ contact and a vacuum-separated slab with three-dimensional periodic boundary conditions. Its reduced cell changes the in-plane coincidence: it is **not** a substitute for the archived 12.51983° low-mismatch interface.

## Recommended candidate: `M77_7site_4Ag_shared_strain`

| Property | Design specification |
|---|---|
| Total atoms | 77 |
| Composition | Ag₂₈Ti₂₁Si₁₄C₁₄ |
| Ti₃SiC₂ slab | Seven atomic planes in the source order Si–Ti–C–Ti–C–Ti–Si; 7 sites per plane = 49 atoms; Si terminated on both faces |
| Ag film | Ag(111), 4 close-packed ABC-stacked layers, 7 atoms per layer = 28 atoms |
| In-plane cell | Hexagonal, index-7 surface repeat; nominal common edge 7.8810 Å |
| Twist | 0° coherent aligned cell; the 12.51983° source twist is not retained |
| Substrate strain | About −3.00% biaxial compression |
| Ag strain | About +3.10% biaxial tension |
| Interface registry | Geometric Si-top starting registry: one bottom-layer Ag site directly above a top-surface Si site per primitive repeat. It is not DFT-ranked. |
| Vertical geometry | Preserve the source's seven Ti₃SiC₂ plane heights and four Ag layer spacings; source Ag layer spacing is about 2.359 Å. Use the source plane-gap reference 2.40 Å. |
| Distance series | `d_ref−10% = 2.16 Å`, `d_ref = 2.40 Å`, `d_ref+10% = 2.64 Å` |
| Vacuum/PBC | Nominal 19.6 Å vacuum above the top Ag layer; PBC = T,T,T, as in the source slab |

Here `d_ref=2.40 Å` is the geometric plane gap from the existing `AGTSC_SI_GAP_2p4A` input. It is **not** an optimized or experimentally established equilibrium distance. The three future scan geometries would differ only by rigid translation of the complete Ag film along z. The substrate, in-plane cell, registry, atom order, and Ag internal stacking would stay fixed across the scan.

The area-7 hexagonal supercell gives a larger lateral repeat than a 2×2 cell while keeping the total below 80 atoms. The shared-strain value is the geometric mean of the inferred unstrained Ag and source Ti₃SiC₂ lattice parameters, which approximately splits the lattice mismatch in logarithmic strain. This is a geometry compromise, not an elasticity-optimized partition. It preserves the layer-resolved coordination topology and keeps the maximum in-plane strain near 3.1%, but changes both materials' in-plane bond lengths and loses the 12.51983° low-mismatch twist. The candidate is only a low-cost local-contact screening model; it cannot establish the energetics of the original low-strain interface or be energy-ranked directly against it.

An alternate `M77_TSC_FIXED` variant retains the Ti₃SiC₂ lattice (0% substrate strain) and puts the full mismatch into Ag (+6.29%). It has the same 77 atoms and an 8.1251 Å edge. This is useful only if preserving the unstrained substrate is prioritized over Ag local bond lengths. Neither strain split is validated; selecting between them needs an explicit modeling decision before coordinates are frozen.

## Compact alternative: `M52_2x2_6Ag`

An optional lower-count alternative is a 2×2 basal cell with four substrate sites per atomic plane (28 Ti₃SiC₂ atoms) and six Ag(111) layers with four atoms per layer (24 Ag), for 52 atoms total. It uses the same Si-top registry and the same nominal `2.16/2.40/2.64 Å` gap series. Its nominal in-plane edge is 6.1420 Å; preserving the substrate lattice requires about +6.29% Ag tension. Six Ag layers preserve a film thickness near the archived four-layer film, but the smaller lateral repeat increases periodic-image coupling. It is a fallback for resource triage, not the preferred local-environment model.

## Structural source and geometry checks required before freezing coordinates

The design derives its seven-plane sequence, basal spacing, Ag layer spacing, gap convention, and vacuum from the archived 217-atom geometry family. The 217-atom input reports 19 Ti₃SiC₂ sites per plane, 21 Ag atoms per layer, 12.51983° twist, +1.1029% Ag strain, 0% substrate strain, and about 20 Å upper vacuum. From the 19:21 area-site ratio and stated Ag strain, the unstrained Ag in-plane spacing is inferred as about 2.8892 Å and the Ti₃SiC₂ spacing as 3.0710 Å. For the recommended area-7 cell, an integer 2D supercell transform with determinant 7 is required; all seven positions in each plane must be reconstructed from the crystallographic parent, not selected by truncating the 217-atom coordinates.

Before any structure is frozen, an offline geometry-generation review must verify:

1. Correct element counts and the Si–Ti–C–Ti–C–Ti–Si substrate stacking.
2. Hexagonal in-plane periodicity and identical cell vectors for the Ag and substrate portions.
3. Ag(111) ABC stacking and four complete Ag layers.
4. The exact top-Si registry and plane-gap definition for all three distances.
5. Minimum periodic Ag–Si, Ag–Ti, Ag–C, and same-phase distances; no unintended overlaps.
6. Vacuum and PBC are retained; no atom/cell relaxation is performed during geometry construction.
7. Every generated coordinate file has a parent/source manifest and SHA-256.

No coordinates are generated in this task, so neither design is yet an executable DFT input. A later resource preflight must measure actual memory; reducing atom count alone does not certify that GPAW will fit within 32 GiB.

## Method boundaries

The existing 15-point input package remains the complete original proposal: six bulk structures and nine 217-atom interface-gap structures. This design does not replace, rescale, or relabel any of them. No result from an unrelaxed minimal geometry may be treated as a DFT label before an authorized, converged SCF and result audit.
