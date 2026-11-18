# Million-atom ACE static check

This is a capability check for the MACE-derived ACE sensitivity model. It reads the repository's 1,021,793-atom QSGS data file and evaluates one LAMMPS `run 0` on four MPI ranks. It does not integrate atoms, heat the system, or measure sustained MD throughput.

## Result

- LAMMPS reached `MACE_ACE_MILLION_STATIC_COMPLETE` in 35 s wall time.
- Potential energy: −6,979,106.6 eV.
- Maximum initial force: 10.336 eV/Å.
- Average neighbor count: 154.7 per atom.
- LAMMPS reported 204–214 MiB memory per MPI rank.

## Force distribution audit

A second `run 0` screened per-atom forces and found 330,452 atoms above 0.5 eV/Å (32.3%), 203,536 above 1 eV/Å (19.9%), and 2,143 above 5 eV/Å (0.21%). The minimum nearest-neighbor distance is 2.088 Å, with no pair below 1.5 Å. High forces are spread across the material-bearing region, rather than arising from an obvious atomic overlap. These are model forces on an unrelaxed QSGS geometry, not fit residuals.

The detailed counts, z bins, and top-force atoms are in `force_audit_summary.json`. Recreate the filtered force dump using `in.mace_ace_force_audit`; the run does not integrate atoms.

These checks establish million-atom loading and static evaluation only. They do not establish that the structure/model is ready for MD.

## Reproduce

```bash
export PATH=/workspace/.venvs/qsgs/bin:$PATH
export LD_LIBRARY_PATH=/workspace/.venvs/qsgs/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mpiexec -np 4 /workspace/potential-comparison/software/lammps-mpich-build/lmp \
  -in in.mace_ace_static -log lammps.log -screen screen.log
```

The `.yace` model is the 0.05-force-weight `Ti_gap` sensitivity candidate. Dataset and model SHA-256 values are recorded in `../ARTIFACT_SHA256SUMS.txt` and `benchmark_summary.json`.
