# AG-SI_INTERFACE_DEV_SOURCE_REVIEW

**Decision: `NO_VALID_SOURCE`.** A credible paper provides three named Ag(111)/Si-terminated Ti₃SiC₂(0001) registry archetypes and scalar equilibrium-spacing/adhesion data. It does not provide a complete machine-readable relaxed structure or the exact Ti₃SiC₂ slab atom count/coordinates needed to freeze a reproducible DEV input. The repository's existing Ag-Si model is a proposal made from a synthetic lateral offset and belongs to an already-scored COD/Kisi family. It must not be promoted to DEV.

No DFT, MACE training, structure optimization, or sampling was performed. Existing labels and their roles remain untouched. The Ag64 data were not used.

## Source review

The strongest reviewed source is Zhang et al., “First-Principles Study of the Interfacial Stability, Electronic Structure and Alloying Effects at the Ti₃SiC₂(0001)/Ag(111) Interface,” *Coatings* 16 (2026) 995, DOI [10.3390/coatings16080995](https://doi.org/10.3390/coatings16080995). The supplied full-text PDF was read locally; its SHA-256 is `fb2a4ad4a659ce7840aa88bb0dcd7f251c02ce7a5fe80441d86309c7d9786ee6`. The PDF itself is not included in this repository report.

The paper establishes:

- The Si termination is one of six explicitly modeled Ti₃SiC₂(0001) terminations (Sec. 3.1; Table 2).
- Three high-symmetry registry archetypes are explicitly defined (Sec. 3.2 and Fig. 4): **OT** (a terminal surface atom directly above an Ag atom), **MT** (bridge midpoint of two adjacent Ag atoms), and **HCP** (threefold hollow at the center of three adjacent Ag atoms). For the Si-terminated case these definitions can guide a deterministic scan; they are not arbitrary offsets.
- The interface model uses a seven-layer Ag(111) slab, a 1×1 lateral cell, Ag coherently strained to the Ti₃SiC₂ in-plane lattice `a=3.071 Å`, and a 15 Å vacuum layer. The paper reports a 6.3% mismatch relative to the Ag(111) nearest-neighbor spacing 2.889 Å.
- Two bottom Ti₃SiC₂ layers are fixed while the remaining atoms, including the Ag overlayer and interface, are relaxed. A 15→20 Å vacuum comparison was included in the paper's convergence discussion, with the representative C(TiC)-HCP adhesion changing by less than 0.1 J/m².
- For the Si-terminated interface, Table 2 reports relaxed `(d₀, Wad)` values: OT `(2.0 Å, 4.38 J/m²)`, MT `(1.8 Å, 5.42 J/m²)`, HCP `(1.7 Å, 6.23 J/m²)`. The paper defines `d₀` as the equilibrium separation from its UBER/relaxation workflow, but does not provide enough coordinate metadata to equate it unambiguously with the project's mean-plane gap convention.

The article and checked supplement do **not** provide a full coordinate file, optimized structure, exact Ti₃SiC₂ substrate slab thickness/atom list, or per-atom force arrays. The text says the substrate layer count was varied to choose an optimum, but the count is not stated in the reviewed methods. Figure 8's four layer-resolved PDOS panels do not prove the total substrate slab thickness. Therefore the article supports registry *categories* and a property-level benchmark, not a ready-to-run geometry or force-labelled DEV parent.

## Existing repository candidate: not a valid DEV source

The available structure `AgSi_cod9009647_slab_v19_dev_proposal.extxyz` has SHA-256 `8fe44bee35e17ebfc433e059925a7a85f0dbb4e775090a75ae1883b3eaa6c36c`. It has 52 atoms (`C16Ag4Si8Ti24`) and comes from COD 9009647/Kisi 1998. The top substrate plane is Si (mean z `23.6058765 Å`), and the single Ag layer is at z `26.2058765 Å`, giving a constructed gap of `2.6 Å`.

Its builder `build_from_cif.py` uses lateral offsets `(i+0.37)a+(j+0.21)b`; this is an implementation choice, not OT/MT/HCP from a source-qualified structure or a relaxed minimum. The same COD-derived 52-atom lineage was used in prior model scoring, and the paper audit found A's three-point slab to be the same atomic geometry after a uniform z translation. This candidate is therefore neither a verified literature registry nor an independent DEV family. It remains a proposal; no labels were changed.

## Minimal deterministic registry and distance scan proposal

Use the published **OT / MT / HCP** archetypes, not random fractional shifts. For each registry use five rigid-gap points centered on that registry's published relaxed `d₀` value, with offsets `−0.20, −0.10, 0, +0.10, +0.20 Å`. Define the project's scan coordinate explicitly as the mean z of the outermost Si termination plane to the mean z of the bottom Ag(111) plane. Because the publication does not map its `d₀` atom-plane convention to coordinates, its tabulated `d₀` values are *initial scan centers only*, not frozen project equilibrium gaps.

This is a 3×5 = **15 fixed-geometry single-point scan design**:

| Registry | Published relaxed d₀ used only as center (Å) | Proposed plane-gap points (Å) |
|---|---:|---|
| OT | 2.0 | 1.8, 1.9, 2.0, 2.1, 2.2 |
| MT | 1.8 | 1.6, 1.7, 1.8, 1.9, 2.0 |
| HCP | 1.7 | 1.5, 1.6, 1.7, 1.8, 1.9 |

Keep each registry and all its distance points in one parent family and one DEV role if later authorized. First compare the scan energies/forces and fit each registry's own gap curve; only then can one select a registry and a project-specific equilibrium distance. Do not call the table's `d₀` an atom-pair distance or use the published Wad values as force labels. The 15-point count is a proposed screening design, not a proven sufficient dataset.

## Candidate geometry settings to carry into a later input-build review

This is a planning freeze only; it does not create a coordinate-level structure release.

| Field | Proposed setting / evidence | Freeze state |
|---|---|---|
| Surface | Ti₃SiC₂(0001), Si terminated; interface orientation Ag(111) | Source-backed choice |
| Registry candidates | OT, MT, HCP as defined in Zhang Fig. 4 | Source-backed categories; exact coordinates absent |
| Ag slab | 7 atomic layers | Source-backed |
| In-plane cell | 1×1 coherent hexagonal cell; `a=b=3.071 Å`, `γ=120°`; strain Ag to substrate | Source-backed paper model; needs project method review |
| Ti₃SiC₂ slab thickness | Obtain exact optimized substrate slab from authors/repository; paper's reviewed methods do not state the final layer count | **Not frozen** |
| z vacuum for project GPAW candidate | 18 Å base with 22 Å comparison, retaining the existing project minimum and +4 Å convergence gate; paper used 15 Å and discussed a 20 Å check | Project proposal only; no convergence evidence for this model |
| PBC / dipole | Project GPAW slab convention `[true,true,false]`, `dipolelayer='xy'`, consistent with existing project slab controls; the paper's 3D periodic supercell implementation is not byte-for-byte reproduced | Project proposal; input-level validation still required |
| Atom numbering | Preserve source substrate order; sort by layer (bottom→top), then fractional `(x,y)`, then species; append Ag atoms by layer (bottom→top), then fractional `(x,y)`; assign immutable consecutive `local_new_id=1..N` and record a species/site/layer map | Numbering rule frozen; actual IDs and `N` **not frozen** until source coordinates and substrate layer count exist |
| DFT labels | For any later approved scan, save native energy, force-consistent free energy, every atom's Cartesian force; save stress diagnostically only | Definition proposal; no calculation authorization |

The 18/22 Å project vacuum proposal is deliberately distinct from reproducing the paper's 15/20 Å CASTEP model. The paper's vacuum convergence result was reported for a representative termination/registry, so it cannot replace a GPAW convergence check for this Si-terminated scan.

## Release gates

`READY_FOR_DFT` can be reconsidered only after a source-qualified full coordinate file or a fully documented, independently reviewed structure reconstruction supplies: the exact Ti₃SiC₂ slab layer count and coordinates; exact OT/MT/HCP fractional translations; seven Ag layers and stacking; cell/PBC and atom/site mapping; source hashes; and a parent-family/geometry overlap audit against the old scored COD-derived Ag-Si family. Then generate and hash all 15 fixed-distance structures and perform the approved k-point, cutoff, and vacuum convergence checks. Until those gates pass, the decision remains **`NO_VALID_SOURCE`** and all DFT remains HOLD.
