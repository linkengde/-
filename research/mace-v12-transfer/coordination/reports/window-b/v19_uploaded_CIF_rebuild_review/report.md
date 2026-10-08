# V19 uploaded-CIF rebuild review

## Decision

The uploaded CIF resolves the **source-content gate provisionally**: it contains an internally consistent COD 9009647 / AMCSD 0013252 attribution and the cited 1998 paper DOI, and ASE expands it to the expected 12-atom conventional Ti3SiC2 cell in space group 194. The remote COD bytes were not fetched, so this review does not establish that the uploaded file is byte-identical to the live COD record.

The three slab files reproduce byte-for-byte from the committed CIF and builder in an isolated temporary checkout. They are geometrically consistent, label-free development proposals. They remain one shared bulk-source family and are not equilibrium surfaces. A fixed-ion numerical sensitivity pilot is reasonable after A freezes the surface registry and registers the calculation owners. This review does not authorize DFT.

## Source and reconstruction evidence

The CIF SHA256 is `6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3`. Its internal fields identify COD 9009647, AMCSD 0013252, Kisi et al., “Structure and crystal chemistry of Ti3SiC2,” *Journal of Physics and Chemistry of Solids* 59 (1998), 1437–1443, DOI `10.1016/S0022-3697(98)00226-1`. The CIF gives P63/mmc (#194), `a=b=3.0575 Å`, `c=17.6235 Å`, and `Z=2`. It has four asymmetric-site rows: Ti(2a) `(0,0,0)`, Ti(4f) `(0.66667,0.33333,0.13550)`, Si(2b) `(0,0,0.25)`, and C(4f) `(0.33333,0.66667,0.07220)`. The expanded sites are C4Si2Ti6 (12 atoms); the multiplicities are 2, 4, 2, and 4. The CIF has no explicit occupancy column; ASE's full-occupancy default is consistent with the stated formula, and no defect occupancy is inferred.

The old prototype values `a=3.07 Å`, `c=17.68 Å`, Ti `z=0.135`, and C `z=0.567` differ from the CIF cell lengths by −0.407% and −0.320%. Comparing the full symmetry-expanded, species-resolved 12-site sets with a single fractional origin shift of `(0,0,0.5)` gives RMS deviation 0.053 Å and maximum matched deviation 0.092 Å in the CIF metric. This treats the C `z` coordinate as an orbit/origin representative; comparing `0.567` directly with `0.07220` would be misleading. The bulk parent is therefore consistent with the old prototype up to the origin convention and small parameter changes.

All seven entries in the source bundle's SHA256 inventory pass. An isolated rerun of `build_from_cif.py` reproduced all three extxyz inputs and `manifest.json` byte-for-byte. Each candidate has 48 CIF-derived substrate atoms plus four Ag atoms, formula C16Ag4Si8Ti24, full PBC, unique IDs in its own 19.1-million range, only a two-atom `central_pair` marker, and no energy, force, or calculator fields. Ordered substrate positions match the CIF expansion, chosen termination cut, 2×2×1 repeat, and +6 Å z translation to below `1e-8 Å` MIC deviation.

## Geometry and changed registry

| Proposal | Substrate top / bottom | Total top layer | Ag gap | Vacuum from atom z extent | Minimum marked Ag–X distance | Minimum framework C–Ti |
|---|---|---|---:|---:|---:|---:|
| AgC | C / Ti | Ag | 2.40 Å | 15.25 Å | 2.526 Å | 2.088 Å |
| AgSi | Si / Ti | Ag | 2.60 Å | 16.00 Å | 2.780 Å | 2.088 Å |
| AgTi | Ti / C | Ag | 2.60 Å | 15.10 Å | 2.980 Å | 2.088 Å |

The manifest's `top_surface_species` means the substrate termination before Ag is placed. In the full structures, Ag is the outermost layer. The one-sided Ag cover, different bottom termination, and one conventional-cell slab thickness make each cell asymmetric; forces and energies can depend on vacuum and slab-image electrostatics. The framework C–Ti minimum is source-derived and should be judged as a C–Ti neighbor, not by an undifferentiated all-pairs cutoff.

Ag in-plane nearest spacing is 3.0575 Å, equal to the substrate repeat `a` by construction, so the imposed Ag/substrate mismatch is 0%. That does **not** mean zero elastic strain relative to free Ag. For scale only, using the commonly quoted fcc-Ag nearest-neighbor spacing of about 2.89 Å would imply roughly +5.8% tensile mismatch; that Ag reference was not independently sourced in this audit, so it is not a release-quality strain value. The earlier B audit JSON also reported 100% by mistakenly comparing the 6.14 Å 2×2 supercell vector with primitive `a=3.07 Å`; the primitive Ag–Ag spacing is the correct comparator.

Although the CIF and old assumed prototype share an origin-equivalent bulk parent, the complete interfaces are not mere relabelings. Recomputed nearest Ag–substrate distances change from the old proposal by −0.284 Å (AgC), −0.001 Å (AgSi), and +0.262 Å (AgTi). The two large changes show that the Ag registry relative to the termination changed when the CIF coordinate convention was used. A should freeze and explicitly approve these registries before any DFT; IDs or the matched bulk cell alone do not establish interface equivalence.

## Isolation against permitted train and development geometries

The geometry-only audit used the 36-frame V19 training candidate and only the three explicitly named V17 development input geometries; reference energy/force columns were not interpreted. The 36 training frames have four composition families: C13Ag2Si10Ti23 (12), C8Ag10Si2Ti12 (1), C13Ag17Si8Ti13 (14), and C11Ag28SiTi9 (9). The reused V17 development frames are 48-atom Ag2C13Si10Ti23, 51-atom Ag17C13Si8Ti13, and 49-atom Ag28C11SiTi9 inputs. None matches the new C16Ag4Si8Ti24 composition. A position-based periodic-contact fingerprint at 3.3 Å (descriptive only, not a chemical gate) found zero short periodic-image pairs in all three V17 development inputs and in 35/36 training frames; the remaining compact training frame has 25. The new slabs have 83–92 such pairs, dominated by framework contacts across the repeated hexagonal plane. Thus the geometry/topology differs from the audited sparse structures; it does not establish broad morphology coverage.

All three proposals share the same CIF bulk parent and composition; they differ by substrate termination and Ag registry. They are therefore related development probes, not three independent source families. Their `CIF_DERIVED_DEVELOPMENT_PROPOSAL` role is appropriate pending A's review; they must not be promoted to blind-test evidence merely because they are new files or have new IDs.

## Numerical pilot recommendation

The AgSi plan in `v19_new_parent_prototypes/convergence_and_source_plan.json` is a reasonable minimal first pilot after A reviews the CIF-derived AgSi geometry and registers the jobs: same frozen input, PW-PBE 500 eV, Fermi smearing 0.1 eV, fixed ions, and matched Gamma 1×1×1 versus 2×2×1. Record native extrapolated energy and free energy separately, full force-vector differences, maximum per-atom force change, and the marked-pair projection. Use the proposed sensitivity budgets only as convergence triggers: 2 meV/atom, 0.01 eV/Å force-vector RMSE, and 0.02 eV/Å pair projection. If exceeded, try 3×3×1; two grids alone do not prove convergence.

At the selected in-plane grid, the smallest vacuum check is the identical AgSi input with atom positions and in-plane cell unchanged and the cell extended only enough to give 20 Å vacuum (about 4.00 Å added to the present c length). Compare the same energy and force observables. Because the slab is asymmetric, check whether the GPAW 26.7 plane-wave setup supports a z-dipole correction under this boundary setup before proposing an on/off comparison; do not guess an API or correction. If supported, make it a matched pair at the 20 Å cell. If not, report that limitation and use vacuum sensitivity as the quantified test.

If Gamma fails the stated sensitivity budgets, do not silently use the denser-grid result as “truth” against Gamma-trained targets. Freeze one target convention for subsequent labels and either keep this pilot as a convergence diagnosis or provide A a common-target/relabel plan before model evaluation or dataset integration. The existing plan's network-blocked source status is now superseded by the uploaded CIF content, but its surface/geometry and numerical gates remain open.

## Scope incident

During the initial exploratory scan, a broad `*.extxyz` glob under `v17_split_geometry_design/final_set/geometries` also reached test-named inputs. A temporary parser split their atom rows and retained species plus xyz to form composition summaries; those geometries were not used in the final comparison or candidate decision, and no energy/force values were evaluated or written to report artifacts. This nevertheless accessed test geometry and violates the task's no-sealed-geometry boundary. The final `audit.py` uses only the three whitelisted V17 development inputs and the training36 file. A should treat this B session as having accessed test geometry and decide whether those V17 test structures remain eligible as sealed evidence.

No DFT, MACE inference/training, MD, or TTM was run. The report, machine-readable audit, and reproducible script are under this directory; `SHA256SUMS.txt` inventories them.
