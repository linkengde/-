# Stage E5 tied-block mapping audit

The E3 map identifies `A`, `B`, `lambda1`, and `lambda2` as equal-valued radial fields shared by each unordered first/second species pair across the third-element triplet slots. E4's mapping audit also checked the matching LAMMPS 22 Jul 2025 Update 4 `pair_tersoff.cpp` implementation (SHA-256 `5edfb3b611e3cc6d07f7cc894eb7354475b607f5e0db87d510a0bd62ceb08d93`): radial pair terms read `[i,j,j]`, while the angular/bond-order paths use triplet entries. Accordingly, this test scales all six ordered triplet rows in one unordered-pair field block together. No individual triplet row was perturbed.

Original TSCBOP SHA-256: `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`. E4 source map SHA-256: `a220e6a05a486653b304ab6c69d32e8a33e78d3effed6acfb1ce1095d8630355`. The original file was copied into `potentials/TSCBOP.source.tersoff`; the source path was only read. E3 reports that paper parameters were externally compared, while E-local access to the paper was not verified; this audit does not assert independent paper verification or establish fit bounds.

| Block | Field | Pair | Base value | Tied table rows | −1% value | +1% value |
|---|---|---|---:|---:|---:|---:|
| A_Ti-C | A | Ti-C | 593.352371 | 6 | 587.41884729 | 599.28589471 |
| B_Ti-C | B | Ti-C | 389.426084 | 6 | 385.53182316 | 393.32034484 |
| lambda1_Ti-C | lambda1 | Ti-C | 1.940953 | 6 | 1.92154347 | 1.96036253 |
| lambda2_Ti-C | lambda2 | Ti-C | 1.62941 | 6 | 1.6131159 | 1.6457041 |
| A_Ti-Si | A | Ti-Si | 3790.763031 | 6 | 3752.85540069 | 3828.67066131 |
| B_Ti-Si | B | Ti-Si | 230.770049 | 6 | 228.46234851 | 233.07774949 |
| A_Si-C | A | Si-C | 46383.947502 | 6 | 45920.108027 | 46847.786977 |
| B_Si-C | B | Si-C | 1752.384773 | 6 | 1734.86092527 | 1769.90862073 |


Each variant SHA and changed-row count is in `parameter_variants.csv`; `parameter_blocks.csv` records the six line/triplet members and values. The 16 variants each start from the same source file and change only the selected field in its six tied rows. All other Tersoff fields remain unchanged.
