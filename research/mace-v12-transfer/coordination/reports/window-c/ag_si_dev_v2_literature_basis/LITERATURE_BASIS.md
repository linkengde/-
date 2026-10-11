# Ag–Si DEV v2 literature basis

## Verified citation and source file

A. H. M. Abdul Wasey, R. Batabyal, J. C. Mahato, B. N. Dev, Y. Kawazoe, and G. P. Das, “First principles electronic structure of coincidence site epitaxial Ag/Si(111) interface,” *physica status solidi (b)* **250**(7), 1313–1319 (2013). DOI: [10.1002/pssb.201248542](https://doi.org/10.1002/pssb.201248542).

The user-provided primary PDF, `First principles electronic structure of coincide.pdf` (7 pages; SHA-256 `d431db3eed82682835869a8e434c0afd1293a74a888c96fef9a09737f49d8630`), was checked against its title-page DOI and text on pp. 1313–1315, including Section 2, Section 3.1, and Figs. 1 and 3. The PDF is not included in this public repository; this record contains only citation facts and a source-file digest.

## Intended use

This is a literature-supported **design-method reference** for a future Ag/Ti₃SiC₂ interface whose MAX-phase surface is Si-terminated and oriented as Ti₃SiC₂(0001), with an Ag(111) overlayer. The cited system is elemental Si(111), not Ti₃SiC₂. Its coincidence cell, registry, and local distances therefore cannot be copied as the target geometry or treated as target DFT labels.

## What the paper supports

1. **Ag(111)/Si(111) coincidence-site epitaxy.** The paper studies an Ag(111) overlayer on Si(111). It reports fcc Ag and diamond-Si lattice parameters of 4.09 Å and 5.43 Å and describes the approximately 25% mismatch and its accommodation by coincidence matching.
2. **4×Ag(111)/3×Si(111) matching idea.** The paper places four Ag(111) surface unit cells on three Si(111) surface unit cells and reports the residual mismatch/effective Ag strain as about 0.3%. This supports considering a coincidence-matching construction, but the in-plane vectors and strain must be recomputed for the actual Ti₃SiC₂ surface cell.
3. **Registry rationale.** In the illustrated model, every fifth Ag atom is placed directly above the fourth Si atom. The paper explains that interface sites have inequivalent subsurface environments and that the Ag layer relaxes/buckles. It presents a physically motivated registry example; it does **not** establish that this is the lowest-energy registry among a systematic scan. It also does not directly prescribe registry for the MAX-phase surface.
4. **Interface-distance reference.** For the relaxed Ag/Si interface, the paper reports Ag–Si bond lengths spanning approximately 2.53–2.72 Å; two inequivalent Ag triangles are reported at approximately 2.69 and 2.72 Å. These are local bond lengths, **not** the mean-plane separation (d_0). They may serve only as qualitative scale/sanity information for designing an independently defined gap scan; do not set the target Ti₃SiC₂ (d_0) equal to these values.
5. **DFT interface-construction method.** The paper describes six Si(111) bilayers, one Ag monolayer, and H passivation of the opposite Si-terminated side; it says the Ag/Si slab coordinates were fully relaxed. A separate figure shows four Ag atomic layers on six Si bilayers. The reported calculations use VASP and PAW potentials for H, Si, and Ag; CA–PZ LDA is used for LDA calculations, with PBE used in some GGA calculations. The methods section reports a Monkhorst–Pack 4×4×2 mesh, a 500 eV plane-wave cutoff, total-energy self-consistency to 0.1 meV, and geometry optimization until forces are below 0.001 eV/Å with conjugate gradients. For a buried Ag/Si supercell, the article separately reports an LDA calculation with a 5×5×1 mesh. These are the paper's settings, not validated settings for GPAW or for the larger Ti₃SiC₂ target interface.

## Limits for the Ag/Ti₃SiC₂ DEV design

The paper supplies a precedent for coincidence matching, one registry example, local Ag–Si bond-length scales, and a slab-relaxation workflow. It does not provide a machine-readable complete target interface geometry, Ti₃SiC₂ slab coordinates, the target interface's (d_0), or complete per-atom energy/force labels. The numerical setup also mixes LDA and PBE calculations and uses VASP, so it is not a drop-in replacement for the project's agreed PBE/GPAW setup.

The actual DEV input still needs a source-qualified Si-terminated Ti₃SiC₂ parent slab, verified layer count and atom map, explicit in-plane matching/strain, candidate registry coordinates, a defined mean-plane gap and scan, vacuum/PBC settings, and structure hashes. The previously reviewed COD-derived 52-atom proposal remains a separate proposal from an already-scored family; this citation neither validates it nor makes it independent.

**Literature status:** `LITERATURE_SUPPORTED_PENDING_STRUCTURE_GENERATION`.

**Execution status:** `HOLD`. No structure was generated, no data role was changed, and no DFT, MD, or MACE training was run; no potential was modified.
