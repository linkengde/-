# Irregular three-dimensional Ag–Ti₃SiC₂ bicontinuous geometry

This is a reproducible approximation of the user's 3-D reference: both phases form interpenetrating, continuous volumes through the model. It uses independently generated ideal FCC Ag and ideal Ti₃SiC₂ parent crystals classified with a seeded smooth random field; it does not reuse coordinates from a melt trajectory and does not implement the referenced four-parameter generator. A single image cannot uniquely specify a 3-D volume.

The output is a geometry candidate only. It has not been evaluated with MACE V17, relaxed, heated, or validated for force-field stability. The V1 candidate with sub-ångström overlaps is retained as a failed attempt and is not used.

- Atom count: 37,711; Ag mass fraction: 59.0766%.
- Box: 110.000 × 98.085 × 100.921 Å; PBC: X free, Y/Z periodic.
- Cross-phase nearest distance after deterministic interface-site construction: 2.200 Å. This is a geometry check only, not proof of stable forces.
- Ag largest-component fraction: 99.93926%; X span: 50.63–105.85 Å; Y/Z winding: True/True.
- Ti₃SiC₂ largest-component fraction: 99.86253%; X span: 50.63–105.89 Å; Y/Z winding: True/True.
- Interface site selection removed 3,976 boundary lattice sites to resolve initial unlike-phase contacts below 2.20 Å. Exact removed sites are listed in `INTERFACE_SITE_REMOVALS.csv`; this changes local interface coordination and requires subsequent static-force review before MD.

Files: `RANDOM_BICONTINUOUS_V2.data` (LAMMPS atomic style), `RANDOM_BICONTINUOUS_V2.extxyz` (OVITO/ASE), `PHASE_ATOM_GROUPS_V2.csv`, `BICONNECTIVE_3D_ATOM_VIEW.png`, `BICONNECTIVE_YZ_SLICES.png`, and `RANDOM_BICONTINUOUS_V2_AUDIT.json`.
