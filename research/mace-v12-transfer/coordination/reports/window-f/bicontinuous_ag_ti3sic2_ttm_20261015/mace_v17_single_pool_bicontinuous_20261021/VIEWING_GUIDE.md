# Compact 3-D Ag–Ti3SiC2 single-source heating pilot

This preliminary pilot was built independently from the ideal Ti3SiC2 parent crystal and FCC Ag lattice, then assigned with a smooth 3-D random field. It does not reuse atom coordinates from a historical melt trajectory.

- `COMPACT_BICONTINUOUS_VIEW.png`: quick perspective projection; gold is Ag, magenta is Ti3SiC2, red circles mark the single localized Ag heating patch.
- `COMPACT_BICONTINUOUS.extxyz`: open in OVITO; X is nonperiodic with 8 A vacuum at both ends, Y/Z are periodic.
- `COMPACT_BICONTINUOUS.data`: LAMMPS atomic-style geometry file.
- `ATOM_GROUPS.csv`: atom IDs, phases, hot-source, sink and interface membership.
- `GEOMETRY_AUDIT.json`: composition, strain, contacts and periodic connectivity.
- `STATIC_RUN0.extxyz` / `STATIC_RUN0_REPORT.json`: initial MACE V17 frozen static result.

The current pilot contains 766 atoms (247 Ag and 519 Ti3SiC2 atoms), 61.01 wt% Ag, and a 17.62 A material thickness. Both phases are geometrically connected and wind through Y/Z periodic boundaries. The hot-source patch contains 41 Ag atoms. The TSC sink atoms remain mobile.

This is a CPU-sized, thin bicontinuous pilot. It does **not** yet include the five physical zones shown in the user's reference diagram (fixed support, thermostat, active reservoir, absorbing layer, and vapor space), and is not a full-size model. Its ongoing trajectory is deliberately not included in this upload. A separate zoned geometry draft is being built to represent that diagram.

The static force and geometric checks are exploratory checks only. They do not validate MACE V17 or establish melting. The 39,486-atom source model and all historical files remain untouched.
