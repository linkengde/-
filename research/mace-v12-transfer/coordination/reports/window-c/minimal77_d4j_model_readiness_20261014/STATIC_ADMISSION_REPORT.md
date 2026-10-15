# Static admission report — M77 Ag/Ti₃SiC₂ candidates

## Decision

The three coordinate candidates pass the geometry-only identity and contact checks listed below. The requested LAMMPS `run 0` force audit was **not run**, so maximum-force admission remains **BLOCKED / NOT MEASURED**. This is not a potential pass.

No DFT, geometry optimization, MD, TTM heating, or parameter fitting was run. The coordinates are fixed, unrelaxed candidates only.

## Geometry check

The candidates each contain 77 atoms with composition Ag₂₈Ti₂₁Si₁₄C₁₄. They retain seven Ti₃SiC₂ atomic planes in bottom-to-top order Si–Ti–C–Ti–C–Ti–Si (7 atoms per plane) and four Ag(111) planes in ABC-A stacking (7 atoms per plane). The Ag bottom A plane is registered over the top Si plane. Atom IDs 1–49 are the substrate, ordered by plane then deterministic in-plane site; IDs 50–77 are Ag, ordered by layer then site. The complete mapping is in `ATOM_ORDER.csv`.

| Gap (Å) | Cell edge (Å) | c (Å) | PBC | Top vacuum (Å) | Minimum all-pair distance (Å) | Minimum Ag–Si (Å) | Minimum Ag–Ti (Å) | Minimum Ag–C (Å) | Pairs <1.5 Å |
|---:|---:|---:|:---:|---:|---:|---:|---:|---:|---:|
| 2.16 | 7.880987 | 37.945160 | TTT | 19.840160 | 2.053774 | 2.160000 | 4.529915 | 5.584839 | 0 |
| 2.40 | 7.880987 | 37.945160 | TTT | 19.600160 | 2.053774 | 2.400000 | 4.752820 | 5.813646 | 0 |
| 2.64 | 7.880987 | 37.945160 | TTT | 19.360160 | 2.053774 | 2.640000 | 4.977315 | 6.043322 | 0 |

The minimum all-pair distance is the Ti–C framework distance, not an Ag contact. The geometrical screen finds no sub-1.5 Å contact. This cutoff only screens obvious overlaps; it does not establish a favorable interface or dynamical stability. The Ag–Si distances are exactly the prescribed plane gaps because the starting registry is Si-top.

The in-plane shared primitive length is 2.978733 Å (index-7 hexagonal repeat). Relative to the design references, the applied biaxial strains are −3.0045% for Ti₃SiC₂ (3.071000 Å reference) and +3.0975% for Ag (2.889238 Å reference). The coherent cell has 0° twist and is not a replacement for the archived 12.51983° interface. The 2.40 Å gap is a geometric reference, not an optimized equilibrium distance.

The seven substrate plane heights are copied from the archived 217-atom fixed-geometry interface parent. In-plane fractional site offsets and stacking are reconstructed from COD 9009647 (Kisi, Crossley, Myhra & Barsoum, *J. Phys. Chem. Solids* 59 (1998) 1437–1443, DOI `10.1016/S0022-3697(98)00226-1`); the shared-strain metric and index-7 cell follow the frozen M77 design. This makes the generation provenance explicit: the coordinates are a derived coherent screening geometry, not an author-provided CIF or relaxed interface.

## LAMMPS force gate

`run 0` and maximum force were **not measured**. The current workspace has no `lmp`/`lammps` executable, no LAMMPS shared library reported by `ldconfig`, and no Python `lammps` module in either active `/workspace/.venvs/mace` or `/workspace/.venvs/gpaw-mpi` environment.

The repository contains the old EAM+Tersoff+Morse calibration-v0 components:

| File | SHA-256 | Audit status |
|---|---|---|
| `research/mace-v12-transfer/coordination/reports/window-b/current_heating_potential_handoff_20261009/potentials/Ag_u3.eam` | `cb3d605cd27ec6efb50c4dd7462c2ab3d967e63775dc4400de9134315a44d807` | Original Ag EAM file present; this by itself is not a four-element interface potential. |
| `research/mace-v12-transfer/coordination/reports/window-b/current_heating_potential_handoff_20261009/potentials/TSCBOP.tersoff` | `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7` | Ti–Si–C-only Tersoff file present. |
| `research/mace-v12-transfer/coordination/reports/window-b/current_heating_potential_handoff_20261009/potential.interface_calibration_v0.inc` | `34deced5e08174ba369290b87f5c3bfd2b284b272991adf302b2b0d28c12815f` | Old Morse cross terms; explicitly labeled “NOT YET A PRODUCTION POTENTIAL.” |

The calibration-v0 README says its Stage69 identity is unconfirmed and the heating branch was paused after abnormal cross-potential energy release. Existing fixed-geometry comparisons also report large force mismatch for the old hybrid. Accordingly, even with a LAMMPS runtime, that file set cannot be represented as a validated D4J production potential. The workspace has no separate, validated D4J four-element potential package or D4J-specific LAMMPS executable to bind to this audit.

The exact maximum force, force components by species, and `run 0` convergence status therefore remain unavailable. No substitute calculator or potential was used.

## Result

- Geometry identity and no-hard-overlap screen: **PASS**.
- LAMMPS `run 0` maximum-force gate: **BLOCKED / NOT MEASURED**.
- Full static admission: **HOLD**.

