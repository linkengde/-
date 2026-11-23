# DFT attempt 002: 64-atom liquid Ag candidate

This is a partial record of the second DFT attempt for the density-corrected candidate. It is **not a usable label**. The run used GPAW 26.7.0, PBE/PW500 eV, 2x2x2 k-points, 0.10 eV Fermi smearing, and four MPI ranks. It was deliberately stopped after only one SCF row had been recorded; the SCF was unconverged, no checkpoint was created, and no forces were extracted. The first row energy was -176.119931 eV/cell and eigst log10 value -0.32.

Stopping was a project-priority decision to proceed with an exploratory potential screen, not a convergence success or a code failure. Files preserve input, preflight, execution script, progress metadata, and the raw log.
