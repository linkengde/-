# Ti3SiC2 pyramidal `(10-11)` geometry candidate

Geometry checks passed. This package contains no DFT labels, no potential scoring, and no data-role assignment.

The source is COD 9009647 (SHA-256 `6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3`). Pymatgen 2026.9.24 generated a symmetric Ti/Ti slab using ASE Miller indices `(1,0,1)`, equivalent to hexagonal `(10-11)`, with minimum slab thickness 15.0 A and requested minimum vacuum 18.0 A. ASE 3.29.0 wrote the extxyz files.

Pymatgen reported normalizing four near-ideal fractional coordinates in the CIF. Direct ASE-versus-Pymatgen parent matching gives a maximum site displacement of 1.765e-05 A and RMS 1.441e-05 A; this is recorded as a parser-level difference, not a new source.

The selected parent contains 69 atoms with composition Ti33Si12C24. Its cell lengths are `[3.0575, 17.886756511452823, 39.7475]` A, PBC is `[true,true,true]`, actual surface-normal span is 14.674850 A (requested minimum was 15.0 A), and periodic vacuum projected along the surface normal is 19.365425 A. This non-bulk-stoichiometric symmetric termination is an as-cleaved neutral-atom slab; no surface energy or thermodynamic stability is claimed.

The response geometry moves the inversion-paired terminal Ti source IDs 6 and 62 outward by 0.02 A along the surface normal. Inversion mapping errors are 7.378e-15 A (reference) and 7.378e-15 A (response). Minimum distances are 2.088203 A (reference) and 2.088203 A (response).

This is a new, not-yet-scored geometry family, but it shares the COD 9009647 crystallographic root with existing datasets. It is only a candidate for facet-held-out development diagnostics, not an independent-source blind test. Both candidates are unassigned; TEST remains outside E's package.
