# Legacy hybrid LAMMPS static run-0 evidence

This directory contains three independently executed static evaluations of the frozen calibration-v0 hybrid potential at the exact 135-atom DFT geometry. `full/in.run0` is the original static input with only a higher-precision thermo print format added. `Ag_only/in.run0` and `TSC_only/in.run0` add `delete_atoms group ... compress no` after `read_data` to remove the opposite phase while preserving original LAMMPS IDs. Each uses `run 0`; there are no minimization commands, timestep integration, MD, or TTM.

The source potential files and `structure.data` are byte copies of the audited source directory. The original repository files were not modified. The LAMMPS wheel required `UCX_TLS=tcp,self` in this restricted container; this changes MPI transport selection only, not the physical model. `units metal` gives energy eV and force eV/Å. Source potential/include hashes are recorded in `LAMMPS_ENVIRONMENT.json` and the parent report.
