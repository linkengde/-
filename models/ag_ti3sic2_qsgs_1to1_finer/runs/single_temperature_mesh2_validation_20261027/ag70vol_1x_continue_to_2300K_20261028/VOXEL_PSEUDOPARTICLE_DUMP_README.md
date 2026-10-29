# Ag70 voxel pseudo-particle dump

`ag70_voxel_pseudoparticles_18p636ps.lammpstrj` is a LAMMPS custom-dump-style
file made from the 18.636 ps 1× enthalpy checkpoint. It contains one record for
each native grid voxel, positioned at the voxel centre. These records are
visualization glyphs for a continuum calculation; they are **not atoms**, and
the file is not a molecular-dynamics trajectory.

The columns are `id type x y z temperature liquid_fraction`. Coordinates and
box lengths use ångströms; temperature uses kelvin. The integer dump step is a
frame label in femtoseconds (`time_ps × 1000`), not an MD integration step.

| Type | Meaning |
| --- | --- |
| 1 | Ti₃SiC₂ scaffold |
| 2 | solid Ag (`liquid_fraction ≤ 1e-12`) |
| 3 | partially liquid Ag in the sphere-render display bin (`0 < liquid_fraction < 0.5`) |
| 4 | Ag in the sphere-render liquid display bin (`liquid_fraction ≥ 0.5`) |

The liquid threshold matches `render_atom_style_phase_3d.py`. The export
includes the full 100 × 100 × 50 domain, including voxels hidden by the image's
cutaway. Boundary labels in the dump are periodic x/y and fixed nonperiodic z.

Recreate the dump from this directory with:

```bash
python export_ag70_voxel_dump.py
```

The adjacent JSON file records phase counts, geometry/checkpoint hashes and
the exported dump's SHA-256 digest.
