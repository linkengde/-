# Ti₃SiC₂ bulk lattice sanity check with repository TSCBOP

This is a small, independent static check of the repository's Ti–Si–C TSCBOP file against a periodic CIF-derived Ti₃SiC₂ crystal. It separates parent-crystal lattice compatibility from the Ag interface and QSGS-created free surfaces. It does **not** validate melting, welding, the Ag interface, or the full composite.

## Result

- Starting crystal: COD 9009647 CIF, converted to an orthorhombic cell with matrix `[[1,0,0],[1,2,0],[0,0,1]]`, then repeated 3×3×2 (432 atoms; `C144Si72Ti216`). Fully periodic in x/y/z.
- With the source-CIF cell fixed, initial maximum force magnitude was 2.67488 eV/Å and pressure was 48,513 bar. A FIRE coordinate minimization converged in 49 iterations to maximum force magnitude 0.000685 eV/Å, but pressure remained 42,930 bar.
- Allowing anisotropic orthogonal box relaxation to zero normal stress, CG converged in 293 iterations. Final maximum force magnitude was 0.00585 eV/Å; pressure was 20.8 bar, with diagonal stress components 31.8, −12.2, and 42.9 bar. The cell expanded by 0.774%, 0.776%, and 0.492% along x/y/z relative to the source-CIF supercell.

This confirms that the TSCBOP file has a relaxed basin for the periodic parent crystal, while its zero-stress lattice is modestly larger than the unrelaxed source-CIF cell. It does **not** explain all large forces in the QSGS composite. The existing composite force audit also finds large forces on Ti/C atoms away from Ag, consistent with the many exposed/cut skeleton surfaces and local coordination changes; cross-interface Morse forces remain an additional source. The bulk result therefore narrows the diagnosis but does not make the QSGS model a stable 300 K starting state.

No source potential or composite structure was changed. No MD, thermostat, TTM, or heat source was used. The zero-pressure endpoint is a diagnostic reference only and must not be substituted into the composite.

## Reproduction

From `models/ag_ti3sic2_qsgs_1to1_finer/`:

```bash
OMP_NUM_THREADS=4 lmp -sf omp -in potential_validation/ti3sic2_bulk_tscbop/in.fixed_cell_minimize -log potential_validation/ti3sic2_bulk_tscbop/log.fixed_cell_minimize
OMP_NUM_THREADS=4 lmp -sf omp -in potential_validation/ti3sic2_bulk_tscbop/in.zero_pressure_box_minimize -log potential_validation/ti3sic2_bulk_tscbop/log.zero_pressure_box_minimize
```

The archived LAMMPS data file is sufficient to run both checks. To regenerate it from the repository CIF, install ASE in a separate Python environment and run:

```bash
python potential_validation/ti3sic2_bulk_tscbop/build_bulk_probe.py
```

Checksums for the inputs, logs, outputs, source CIF, and unchanged potential are in `SHA256SUMS.txt`.
