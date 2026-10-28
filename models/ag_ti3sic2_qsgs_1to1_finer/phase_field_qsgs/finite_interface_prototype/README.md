# Finite-interface prototype

This directory contains an isolated conservative finite-volume conduction kernel and its small numerical validation suite. It does not alter the production enthalpy solver or any saved geometry/checkpoint.

Run the tests from the repository root:

```bash
python -m unittest -v models.ag_ti3sic2_qsgs_1to1_finer.phase_field_qsgs.finite_interface_prototype.test_finite_interface_conduction
```

The kernel accepts a 3D phase-label array, physical cell pitches, cell conductivity and volumetric-capacity fields, and a finite `G_K` in W/(m²·K). Fixed-temperature outer faces and periodic axes are optional. Unlike-phase links use the series half-cell plus interface resistance; pairwise heat exchange is equal and opposite.

`G_K=0`, finite `G_K`, and `G_K=∞` are tested against a one-dimensional Ag/Ti₃SiC₂ bilayer series-resistance solution. The transient test uses the current 300 K candidate properties and closes the insulated-domain energy ledger. Test `G_K` values are numerical fixtures only, not interface-property estimates.

See:

- `FINITE_INTERFACE_MODEL_PREPARATION.md` for equations, tests, and limits.
- `MATERIAL_SOURCE_AUDIT.md` for candidate-property provenance and unresolved records.
- `CONTACT_HEAT_SOURCE_AUDIT.md` for missing experimental source inputs.
- `GEOMETRY_AND_RVE_AUDIT.md` for the existing 70:30 raster statistics and micro/nano RVE plan.

