# Exploratory GPAW single point

`run_01/` stores one fixed-geometry GPAW PBE calculation on the 12-atom reconstructed Bai2010 CIF. The run uses 500 eV, Fermi smearing 0.10 eV, 12×12×2 k points, and four MPI ranks. This k-point mesh approximately matches the reciprocal spacing of the existing 2×2×1 Ti₃SiC₂ numerical control at 6×6×2; it does not establish convergence.

The result is diagnostic only. Native and force-consistent free energies are kept separately. The SCF log is losslessly gzip-compressed after the run; `gpaw.log.raw.sha256` records the original uncompressed log hash. The structure and labels are excluded from all train/dev/test datasets, and no relaxation is performed. Re-run only from a fresh checkout after checking that `run_01/` is absent:

```bash
GPAW_MPI_BACKEND=cgpaw GPAW_MPI=/lib/x86_64-linux-gnu/libmpi.so.40 \
LD_LIBRARY_PATH=/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu:/workspace/.venvs/gpaw-mpi/lib \
OMP_NUM_THREADS=1 /usr/bin/mpirun -n 4 \
/workspace/.venvs/gpaw-mpi/bin/python run_exploratory_single_point.py
```
