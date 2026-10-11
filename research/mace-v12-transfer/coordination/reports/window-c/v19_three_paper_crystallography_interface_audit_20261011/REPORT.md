# Ti₃SiC₂ crystallography and Ag/Ti₃SiC₂ interface source audit

Audit date: 2026-10-11
Project: Ag–Ti₃SiC₂ unified four-element MACE potential
Source baseline: `b55c576dd83a248e69caee8c6fbb44a4d7cccb59`
Scope: offline source, crystallography, and data-lineage review only.

## Decision

The package supplies sound experimental crystallographic references for the ideal Ti₃SiC₂ prototype and a useful literature benchmark for Ag/Ti₃SiC₂ adhesion. It does **not** supply an independent Ti₃SiC₂ validation parent or complete Ag/Ti₃SiC₂ DFT training data.

The 1967, 1998 COD/Kisi, and 2000 structures are the same `P6₃/mmc` Ti₃SiC₂ prototype. Rawn 2000 is a very close geometric match to the already scored COD/Kisi parent. Jeitschko–Nowotny 1967 has a larger C-coordinate difference but the same Wyckoff sites and layer sequence, and is not a new independent parent for ML splitting. The Zhang 2026 paper reports 18 scalar interface-spacing/adhesion results, but no complete interface coordinates, relaxed coordinate files, or per-atom forces. Treat it as an external property-level literature reference only.

The Ag–Si three-point series from A also belongs to the old scored 52-atom COD-derived family: the A parent has the exact same 52 atomic coordinates and ordering as the old parent after a uniform `+5.0 Å` z translation; its cell is 10 Å longer and z-PBC is off. A’s minus/plus structures rigidly move only the four-atom Ag layer by `−0.20/+0.20 Å`. This is strong parent-family overlap. Neither those three points nor the Rawn/Jeitschko reconstructions qualify as independent DEV/TEST structures.

Accordingly:

- **Usable physical references:** room-temperature Ti₃SiC₂ lattice/site parameters; Zhang’s ideal-interface `d₀` and `Wad` table as a method-qualified literature benchmark.
- **Possible future training candidates:** reconstructed Rawn or Jeitschko geometries only as same-prototype-family geometry proposals, after explicit data-role approval and new, method-consistent DFT labeling. The existing COD numerical-control labels retain their current role. The papers themselves contain no new DFT energy/force labels.
- **Independent DEV/TEST qualification:** none of the new sources meets it. The old COD parent and A’s displacements are already in the scored lineage. Zhang provides insufficient structural data to establish a distinct parent.
- **Current V19 gates:** five historical backbone controls remain numerical controls; the seven proposed new DFT labels and R30/F30/R30+S training remain HOLD. D/E/F’s six liquid-Ag pilots are unchanged; A/B are unchanged; G/H remain on standby. No DFT, MD, or training was run in this audit.

## Source package verification

All individual input files were received as the listed files, rather than as a single archive image. `sha256sum -c SHA256SUMS.txt` passed all seven listed entries. Each of the four assets enumerated in `MANIFEST.json` matches its manifest size and SHA-256. The package checksum-list file itself has SHA-256 `02c52d404b9913f7e3fe22121e0c5b5a50171632a8e63a8402abac25b0f4e92d`; it is not self-listed. The COD input has the same SHA-256 as the repository’s COD 9009647 CIF. PDF page counts were checked with `pdfinfo`; titles, key tables, methods, and data-availability statements were read from the supplied full text.

The three PDF hashes are recorded in `source_file_audit.csv`. The two historical papers remain private audit inputs and were **not** copied to this repository. The 2026 MDPI article identifies itself as CC BY, and the extracted Table 2 is attributed below; its separate data-availability statement says some or all generated/used data, models, or code may be proprietary/confidential and provided only with restrictions. Therefore the article’s license does not establish access to its missing coordinate or force files.

## Paper-derived Ti₃SiC₂ CIF reconstruction

Two new CIFs are full 12-atom conventional hexagonal cells (`Ti₆Si₂C₄`) derived from the paper-reported lattice and asymmetric-site coordinates. They are not author-supplied CIFs. For unambiguous full-cell atom counts, the generated files use P1 with the 12 symmetry-expanded atoms; an audit field records the parent prototype `P6₃/mmc` (No. 194), the original Wyckoff-site assignment, source citation, and occupancy assumption. Coordinates were expanded with the 24 symmetry operations in the COD 9009647 file.

| Source | Experimental details and values used | Reconstruction caveat |
|---|---|---|
| Jeitschko & Nowotny (1967), single-crystal XRD, DOI `10.1007/BF00899949` | `a=3.068 Å`, `c=17.669 Å`; Ti(I) 4f `z=0.1350(15)`, Ti(II) 2a `(0,0,0)`, Si 2b `(0,0,1/4)`, C 4f `z=0.5675(8)`. Difference-Fourier C position is also reported as `z=0.570(4)`. | The paper discusses possible 2a Ti under-occupancy/Ti–Si exchange and C occupancy uncertainty. The CIF is explicitly an ideal stoichiometric (`occupancy=1`) reconstruction, not a claim that the refined experimental occupancies were all exactly full. |
| Rawn et al. (2000), neutron powder diffraction/Rietveld, DOI `10.4028/www.scientific.net/MSF.321-324.889` | `a=3.0654(1) Å`, `c=17.635(1) Å`; Ti(I) 4f `z=0.1348(2)`, Ti(II) 2a `(0,0,0)`, Si 2b `(0,0,1/4)`, C 4f `z=0.57244(8)`. Room-temperature refinement after thermal aging reports Si occupancy `1.00(2)`. | This is the post-heat-treatment room-temperature structure, not a high-temperature configuration. Rawn also reports temperature-dependent apparent Si occupancy in in-situ data and explains its correlation with displacement parameters/preferred orientation; do not interpret this reconstruction as a liquid/high-temperature label. |

Both generated CIFs contain exactly 12 atoms with nominal formula Ti₆Si₂C₄ and preserve the source’s hexagonal setting and Wyckoff-site tags. Fractional coordinates are serialized to 8 decimal places; the resulting Cartesian rounding is below `1×10⁻⁷ Å`. For a true symmetry-expanded structure, the P1 setting avoids an importer applying the P6₃/mmc operations a second time.

## Periodic geometry and local coordination comparison

The CIF and each paper-derived asymmetric site set were expanded to the 12-atom cell using the same 24 operations. Species-resolved periodic atom mapping searched the global origin translation; minimum-image distances were calculated in the pairwise-averaged hexagonal metric to separate internal-coordinate differences from lattice-constant differences. The minimum global origin shift is `(0,0,0.5)` for each comparison. Per-atom mappings and displacements are in `atom_mapping.csv`; aggregate values are in `periodic_structure_comparison.csv`.

| Paper structure vs COD 9009647 | Ti₂a RMS/max (Å) | Ti₄f RMS/max (Å) | Si₂b RMS/max (Å) | C₄f RMS/max (Å) | all-atom periodic RMS/max (Å) | Δa / Δc / ΔV vs COD (%) |
|---|---:|---:|---:|---:|---:|---:|
| Jeitschko–Nowotny 1967 | 0 / 0 | 0.00882 / 0.00882 | 0 / 0 | 0.08294 / 0.08294 | 0.04815 / 0.08294 | +0.343 / +0.258 / +0.948 |
| Rawn et al. 2000 | 0 / 0 | 0.01234 / 0.01234 | 0 / 0 | 0.00423 / 0.00423 | 0.00753 / 0.01234 | +0.258 / +0.065 / +0.583 |

The site-resolved local coordination table uses explicit first-shell cutoffs (`Ti–C < 2.50 Å`, `Ti–Si < 3.10 Å`) that lie between the six-or-three near neighbors and the next shell for these ideal structures. Results are consistent across all three sources: each Ti(II) 2a has six C neighbors, each Ti(I) 4f has three C neighbors, each C has six Ti neighbors, each Ti(I) 4f has three Si neighbors, and each Si has six Ti(I) neighbors; Ti(II) has no Si within 3.10 Å. The full per-site distances and counts are in `local_coordination_comparison.csv`. The different experimental internal coordinates shift near-neighbor lengths by hundredths of an ångström but do not change the coordination topology.

These are same-prototype geometry comparisons, not evidence that three distinct papers produce three independent statistical samples. Rawn’s references cite both Jeitschko–Nowotny (1967) and Kisi et al. (1998); COD 9009647 is the Kisi et al. crystallographic source. Rawn’s cell and internal coordinates are extremely close to COD. Jeitschko’s C z coordinate differs more, but space group, Wyckoff sites, composition, and layer topology still match. All three therefore remain in one conservative Ti₃SiC₂ parent-family group for any ML split.

## Zhang et al. (2026): what the interface paper does and does not provide

Zhang et al., *Coatings* 16 (2026) 995, DOI `10.3390/coatings16080995`, studied six Ti₃SiC₂(0001) terminations and three lateral registries (OT, MT, HCP), reporting 18 unrelaxed UBER-fit and relaxed `d₀/Wad` pairs in Table 2 (PDF page 9). The values, in ångströms and J m⁻², are transcribed in `zhang_table2_interface_benchmark.csv` with attribution. Across the table, unrelaxed `Wad` ranges from 1.82 to 7.18 J m⁻² and relaxed `Wad` from 2.87 to 9.25 J m⁻²; C(TiC)-HCP has the highest relaxed value, 9.25 J m⁻² at `d₀=1.2 Å`.

The method disclosed is CASTEP/Materials Studio 2020 with PBE and on-the-fly ultrasoft pseudopotentials, Gaussian smearing 0.1 eV, 500 eV cutoff, 9×9×1 surface/interface and 9×9×9 bulk meshes, 15 Å vacuum plus dipole correction, and a coherent 1×1 in-plane registry. The Ag surface is laterally strained to the Ti₃SiC₂ in-plane parameter (3.071 Å; stated mismatch about 6.3% relative to Ag surface nearest-neighbor spacing 2.889 Å). Two bottom Ti₃SiC₂ layers are fixed; remaining atoms are relaxed. Reported optimization criteria include energy change per atom below 10⁻⁵ eV, maximum force below 0.02 eV/Å, maximum displacement below 5×10⁻⁴ Å, and residual stress below 0.05 GPa.

The article supplies schematic figures and scalar `d₀/Wad` results, but the checked article and supplement do not provide a complete atom-coordinate listing/CIF/structure archive, optimized coordinates, total-energy outputs, or per-atom force arrays. The supplement is PDOS figures, not structures. A convergence threshold for maximum force is not a force label. Consequently these data cannot be treated as MACE energy/force training labels, and the table’s `d₀` must not be converted to a universal Ag–Ti/Si/C pair cutoff. `d₀` is the authors’ interface-level equilibrium separation measure; without the full termination-specific coordinates its atom-pair mapping is unavailable.

This remains useful as a property-level trend reference for ideal coherent interfaces, subject to the different electronic-structure implementation/pseudopotential and cell model. It is not a direct apples-to-apples comparison with the existing 52-atom GPAW slab: the current slab has four Ag atoms in a single layer on a COD-derived 2×2 framework, whereas Zhang describes a seven-layer Ag slab in a 1×1 lateral registry and termination-dependent Ti₃SiC₂ slabs. No coordinate-level parent similarity or independence can be established for Zhang’s structures from the supplied files.

## Ag–Si family overlap and leakage assessment

The current A DFT archive contains a 52-atom `C16Ag4Si8Ti24` parent and fixed-geometry Ag-layer z shifts at −0.20 and +0.20 Å. Each output reports SCF convergence, verifier PASS, PBE/PW500, 6×6×1 k mesh, 0.10 eV smearing, and native and force-consistent energies separately; these are usable archived DFT labels only in their already assigned role. They do not form a fresh held-out parent family.

An atom-by-atom comparison of A’s parent input against the older `AgSi_cod9009647_slab_v19_dev_proposal.extxyz` gives a maximum residual of `0.0 Å` after subtracting a uniform `(0,0,+5.0 Å)` translation. Composition, ordered symbols/IDs, and the full 52-atom geometry are identical; the A cell adds 10 Å along z and changes z-PBC from true to false. The A ±0.20 Å geometries leave the framework unchanged and displace all four Ag atoms together. The old COD-derived 52-atom slab was already used in V18/base-model force scoring documented in the V18 audit. This creates direct parent-family overlap and a potential train/score leakage path if any members are split separately.

Zhang’s ideal interface class shares the same Ag/Ti₃SiC₂ chemistry and (0001)/(111) orientation, but its described 1×1/7-layer model is not the same geometry as A’s 2×2/single-layer slab. Since coordinates are not supplied, it is impossible to prove a geometrical parent collision—or to establish that it is a genuinely independent DEV source. Classify the paper as a separate scalar literature benchmark, not as a verified independent structural parent. Do not split A’s parent and the old 52-atom parent across train/dev/test; group them together as one COD-derived family.

## V19 data-role judgment and next evidence needed

| Item | Physical reference? | Training-candidate status | Independent DEV/TEST? | Current role/action |
|---|---|---|---|---|
| Jeitschko–Nowotny 1967 reconstructed Ti₃SiC₂ | Yes, historical room-temperature prototype and approximate site coordinates | At most a same-family geometry proposal; no DFT E/F labels, and occupancy uncertainty reduces priority | No | Reference only; keep out of training until separately authorized and labeled |
| Rawn 2000 reconstructed Ti₃SiC₂ | Yes, more precise neutron/Rietveld room-temperature coordinates and lattice | At most a same-family geometry proposal; likely redundant with COD; no DFT E/F labels | No; near duplicate to COD | Reference only; do not use as independent DEV |
| COD/Kisi 1998 Ti₃SiC₂ and Ag–Si parent | Yes | Existing numerical-control/scored roles stay unchanged | No; previously scored family | Do not silently convert numerical controls to training or test data |
| Zhang 2026 Table 2 | Yes, qualified scalar literature reference for ideal coherent interface trends | No training label without coordinates and energy/force data | Not established | Store cited table only; request author/supplementary coordinates if exact validation is needed |
| A parent, −0.20, +0.20 Å Ag–Si labels | Yes, converged current GPAW references for their exact fixed geometries | Role changes require explicit approval | No; same already scored COD family | Retain A’s role; family-group together |

No fresh independent Ti₃SiC₂ DEV parent or Ag–Si DEV parent is established by this package. To create one, obtain a structure whose experimental/source lineage and geometry are demonstrably distinct from the scored COD-derived prototype family (for Ti₃SiC₂), or a fully specified interface model with atom positions, termination, registry, composition, cell/PBC and an independent source parent (for Ag–Si). Any future DFT labels must include converged energies and complete forces under one declared energy/free-energy convention; independent status must be frozen before model scoring.

## Reproducibility and output inventory

- `Ti3SiC2_jeitschko_nowotny_1967_reconstructed.cif` — full 12-atom P1 cell expanded from paper-reported P6₃/mmc sites.
- `Ti3SiC2_rawn_2000_reconstructed.cif` — same format for the Rawn neutron refinement.
- `periodic_structure_comparison.csv` and `atom_mapping.csv` — lattice-aware periodic mapping and coordinate residuals.
- `local_coordination_comparison.csv` — first-shell coordination and neighbor distances.
- `zhang_table2_interface_benchmark.csv` — 18 Table 2 rows, attributed to Zhang et al. (2026), Table 2.
- `ag_si_lineage_audit.csv` — A parent/displacements and old scored COD-derived slab identity.
- `source_file_audit.csv` — private package byte lengths, original hashes, page counts, and DOI identity.
- `SHA256SUMS.txt` — SHA-256 of every committed output other than the checksum file itself.

No uploaded historical PDF, package ZIP, sealed structure/label, raw A/B checkpoint, frozen model, or original training data is included in this directory. No DFT, MD, MACE inference, or training was launched for this audit.
