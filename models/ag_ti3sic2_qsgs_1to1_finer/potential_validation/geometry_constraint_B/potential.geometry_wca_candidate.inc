# EXPLORATORY GEOMETRIC-CONFINEMENT MODEL B — NOT A PHYSICAL Ag/MAX INTERFACE POTENTIAL
# Type 1 Ti, 2 Si, 3 C, 4 Ag. Ceramic is intended to be held fixed in the test inputs.
# Ag-ag remains the repository EAM; Ti-Si-C remains repository TSCBOP.
# Ag-skeleton is WCA-only repulsion. rc reuses calibration-v0 Morse r0 solely as
# a provisional contact distance; epsilon=0.5 eV is an exploration parameter.
pair_style hybrid eam tersoff lj/cut 2.40407837
pair_modify pair lj/cut shift yes
pair_coeff 4 4 eam potentials/Ag_u3.eam
pair_coeff * * tersoff potentials/TSCBOP.tersoff Ti Si C NULL
# WCA: rc=2^(1/6)*sigma; shifted LJ is zero-energy, zero-force at rc.
# No attraction, adhesion or interfacial heat-transfer physics is represented.
pair_coeff 1 4 lj/cut 0.5 2.08210927 2.33708864
pair_coeff 2 4 lj/cut 0.5 2.14179034 2.40407837
pair_coeff 3 4 lj/cut 0.5 1.91771227 2.15255924
