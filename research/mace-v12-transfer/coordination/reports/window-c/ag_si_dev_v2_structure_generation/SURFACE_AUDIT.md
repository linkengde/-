# Ti₃SiC₂(0001) Si-termination surface audit

## Source and lineage

The substrate is reconstructed from the repository's COD 9009647 CIF, attributed there to Kisi et al., *Journal of Physics and Chemistry of Solids* **59** (1998), 1437–1443, DOI `10.1016/S0022-3697(98)00226-1`. Source path:

`research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif`

Source SHA-256: `6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3`.

ASE 3.29.0 expands the CIF's four independent crystallographic sites and symmetry to a 12-atom conventional hexagonal cell (two Ti₃SiC₂ formula units), with `a=b=3.0575 Å`, `c=17.6235 Å`, `α=β=90°`, `γ=120°`. The cell is P6₃/mmc, space group No. 194. The fractional layer sequence in one c cell is Ti–C–Ti–Si–Ti–C–Ti–C–Ti–Si–Ti–C after wrapping at the cell boundary.

**Lineage restriction:** this is the COD/Kisi crystal parent used by prior scored Ag–Si structures. The reconstructed interface geometries are not an independent DEV source family. No TRAIN/DEV/TEST role is assigned or changed.

## Selected slab construction

The slab is a deterministic, unrelaxed cleavage/reconstruction from two successive Si-to-Si Ti₃C₂ blocks of the COD parent. It contains two outer Si termination planes plus one internal Si plane. Both outer substrate faces are Si terminated; the upper face contacts Ag, while the lower face remains a free Si surface. No hydrogen or other species is added.

| Property | Frozen candidate geometry |
|---|---|
| Selected termination | Si on both outer Ti₃SiC₂(0001) faces; Ag contacts the upper Si face |
| Substrate layer order, bottom to top | Si–Ti–C–Ti–C–Ti–Si–Ti–C–Ti–C–Ti–Si |
| Substrate atomic planes | 13 (each has one atom in the 1×1 basal cell) |
| Substrate atom counts | Ti₆Si₃C₄; 13 atoms |
| Substrate thickness | 17.6235 Å between outer Si planes |
| Coherent in-plane cell | `a=b=3.071 Å`, `γ=120°`, from the reviewed Ti₃SiC₂/Ag interface paper's 1×1 model; COD fractional coordinates are mapped into this cell |
| Full XYZ cell | `a=b=3.071 Å`, `γ=120°`, `c=71.268266 Å`; `α=β=90°` |
| PBC | x/y periodic, z nonperiodic (`T T F`) |
| Vacuum margins | 18.0 Å below the bottom Si plane and 18.0 Å above the top Ag plane |
| Interface plane gap | 2.600 Å from mean upper Si-plane z to mean bottom Ag-plane z; a common unrelaxed starting value, not a measured equilibrium `d₀` |

The COD experimental in-plane parameter (3.0575 Å) is scaled by about +0.44% to the 3.071 Å coherent-cell parameter used in the reviewed direct interface study. This is a documented geometry construction choice, not a newly optimized Ti₃SiC₂ lattice.

The asymmetric one-sided interface leaves a free Si-terminated lower face. A later DFT input review must account for that surface, slab-thickness effects, and the project's slab electrostatics; this geometry package does not certify those settings or authorize DFT.

## Ag count and composition

Each candidate has 7 Ag(111) atomic layers, with one Ag atom per common 1×1 cell per layer (full 1×1 coverage). Total composition is Ag₇Ti₆Si₃C₄, 20 atoms. Atom IDs are consecutive and stored in each extended XYZ together with subsystem, layer index, source CIF site index, and c-repeat index.

`READY_FOR_DFT_SCREENING` means that complete candidate coordinates and cell/PBC metadata exist and passed the stated geometric overlap check. It does not mean a DFT calculation is authorized, the structure is relaxed, the interface is stable, or that these structures are independent DEV data.
