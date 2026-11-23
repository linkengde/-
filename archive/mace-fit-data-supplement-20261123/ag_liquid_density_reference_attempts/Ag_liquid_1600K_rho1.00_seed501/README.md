# Ag64 liquid-density candidate: preflight attempt

This folder records a GPAW initialization-only attempt for one 64-atom Ag geometry at the cited 1600 K reference density. The geometry is only a candidate proposal from an EAM trajectory; it has no DFT energy or forces.

The dry-run used GPAW 26.7.0, PBE/PW500 eV, 2x2x2 k points, 0.10 eV Fermi smearing, and four MPI ranks. It remained at `Converting LCAO to pw mode` for 20 min 54 s and was stopped before SCF. No label or checkpoint was produced. The input, configuration, preflight log, scripts, and termination/resource record are preserved here for audit.

Do not treat this attempt as DFT data or include the structure in a fitted dataset. The 2x2x2 protocol needs a more targeted startup-cost diagnosis before another production attempt.
