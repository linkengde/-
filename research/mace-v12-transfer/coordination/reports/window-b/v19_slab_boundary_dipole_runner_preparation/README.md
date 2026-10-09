# Disabled slab boundary/dipole draft

Two exact same-ion/same-cell controls use original16.00026725 A vacuum, pbc=(T,T,F), provisional6x6x1/PBE/PW500/Fermi0.10. Off uses poissonsolver={}; on uses {'dipolelayer':'xy'}. Inputs preserve source decimal bytes, ordered elements, historical lammps_id and central_pair; only source PBC header changes. The periodic baseline is not repeated. No relaxation, hidden translation, width or target change.

Source API evidence in installed GPAW26.7: old/pw/hamiltonian.py passes dipolelayer to DipoleCorrection; dipole_correction.py rejects periodic normal and nonorthogonal axes. New/pw/builder.py forwards Poisson parameters; new/pw/poisson.py requires two periodic directions and uses the nonperiodic axis. For this z-normal cell xy is correct. Native energy (force_consistent=False) and free_energy are preserved separately; changing the returned energy key does not alter forces. Summary method includes actual PBC/Poisson, GPAW26.7, PAW-data1.2.1, PBE setup hashes, cutoff/mesh/width, SCF/mixer/parallel and energy keys. Output serializes actual TTF PBC, never hardcoded TTT.

entry.py rejects disabled launch before importing GPAW/creating output. It enforces source/input/inventory, identities/geometry/PBC/method/Poisson, exact registered owner/label, CPU4 and disk3GiB, package/setup pins. Queue uses a local exclusive lock and active-run/archive guards, four MPI ranks and one BLAS thread; runner enforces MPI/ScaLAPACK, no output overwrite, disk1.5GiB per SCF step. Compact archive excludes state.gpw. Complete inventory, convergence/log/progress, exact identity/PBC/Poisson/energy and finite forces are independently verified.

All launch flags remain false and owner null. A must review/move this pack into a production entry, register exact owners/labels, regenerate inventory and integrate compact Git publication/completion. The draft queue publishes start intent but deliberately returns a verified archive to A's publication workflow; it cannot silently authorize itself. Do not execute synthetic fixture results as references. No DFT was run.

Static validation uses temporary, explicitly NONPHYSICAL energies/forces/logs. It checks disabled direct/queue launch and rejects wrong PBC/Poisson/input hash/IDs/method/source/convergence/convention/nonfinite forces and extraneous inventories. Fixtures do not establish runtime SCF or physical readiness. API evidence is source inspection; actual PW slab execution awaits A registration.

From this directory with installed GPAW venv:

```bash
/workspace/.venvs/gpaw-mpi/bin/python test_rejections.py
```

Generated validation.json changes require regenerating SHA256SUMS.txt afterward. This inventory guards reviewed entry scripts/inputs; runtime calculations/ and __pycache__/ are excluded. Result inventory independently covers the complete archive.
