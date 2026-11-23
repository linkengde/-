# Full-size QSGS exploratory geometry and static EFS audit

Date: 2026-11-23  
Decision: **not ready for 300 K MD, TTM, or melt/quench dynamics**

## Scope

This was an exploratory attempt to retain the original 482,471 Ag/W atoms and QSGS phase assignment while reducing the 0.08 Å Ag–W overlap. The original `ag_w_qsgs.data.gz` was read-only and its SHA-256 remains `204961ff92491b02d3517a172377b997b378f47e74c47557ca9b6e76d25f0`.

The geometry-only preconditioner projected Ag–W, W–W, and Ag–Ag pairs toward 2.54, 2.60, and 2.48 Å, respectively, using simultaneous equal-and-opposite corrections. These are screening targets taken from extrema in the v3 local slab scans; they are not fitted force-field parameters. All atoms and atom types were retained. The candidate mass ratio is 0.999766 (relative error 0.0234%).

## Geometry projection result

After 500 projection rounds, the candidate retained both phases' spanning paths under the 95% connectivity screen. The minimum distances were 2.537 Å (Ag–W), 2.596 Å (W–W), and 2.478 Å (Ag–Ag). However, the algorithm had redistributed contacts into same-phase neighborhoods and did not remove all target-distance violations. This output is retained as `ag_w_qsgs_v3_exploratory_geometry.data.gz`; it is not an accepted structure.

## Spring-free Dai EFS static relaxation

The candidate was minimized in LAMMPS with the Dai 2006 EFS `eam/fs` table, `p p f` boundaries, a fixed 0–20 Å bottom layer, and all remaining atoms mobile. The conjugate-gradient minimizer ran 10,000 iterations and stopped at its iteration cap, not at force convergence.

| Check | Result | Practical screen | Outcome |
|---|---:|---:|---|
| Total atoms / mass ratio | 482,471 / 0.999766 | retain all / error ≤0.1% | pass |
| Largest connected share, W / Ag | 98.24% / 99.68% | each ≥95%, x/y winding, spans z | pass at this cutoff |
| Mobile RMS force | 0.00411 eV/Å | ≤0.02 eV/Å | pass |
| Mobile maximum force | 0.650 eV/Å | ≤0.1 eV/Å | fail |
| Minimum Ag–Ag distance | 2.099 Å | ≥85% pure EFS nearest neighbor (2.458 Å) | fail |
| Minimum W–W distance | 2.412 Å | ≥85% pure EFS nearest neighbor (2.326 Å) | pass |
| Minimum Ag–W distance | 2.267 Å | ≥90% v3 relaxed local-slab minimum (2.290 Å) | fail |
| Atoms above z=200 Å | 1,604 (maximum z 202.387 Å) | none before heating | fail |
| Maximum displacement from projected input | 4.652 Å | inspect for major local rearrangement | fail / severe rearrangement |

The low RMS force is dominated by most atoms being well relaxed; it does not cancel the maximum-force outlier or the collapsed Ag–Ag and Ag–W contacts. The potential minimization also displaced atoms into the nominal vacuum before any heating. Therefore this structure fails the proposed practical MD-start gate even though atom count, mass ratio, and percolation are largely retained.

## Decision and limits

- No 300 K equilibration, TTM, laser pulse, melting, or solidification trajectory was started.
- DFT interface-energy mismatch is not the reason for rejection here. The full-size structure itself remains mechanically and geometrically unsuitable.
- This outcome does not invalidate qualitative graduate-level use of the Dai EFS model once a viable structure and temperature protocol are available.
- The next geometry route should reconstruct atomic neighborhoods locally (or build the phases with a controlled interface) rather than continue this all-pairs projection. Preserve the original QSGS morphology and atom-count target when designing that route.

## Reproduction files

- `project_full_qsgs_v3.py`: geometry-only projection.
- `geometry_projection_report.json`: projection outcome and source/candidate hashes.
- `static_efs_minimize.in`: LAMMPS static-relaxation input.
- `static_efs_minimize.log`: complete 10,000-step LAMMPS log.
- `static_efs_audit.json`: parsed forces, distances, topology, and pass/fail values.
- `analyze_static_candidate.py`: reproduces the machine-readable static audit from the minimized dump.
- `minimized.data.gz` and `minimized.dump.gz`: minimized coordinates and forces.

The uncompressed working files are omitted from the pack; decompress the `.gz` files when needed. The serial partial log is kept separately to document the switch to the four-thread OpenMP run.
