# EXPLORATORY B MODEL FOR A RIGID GEOMETRIC SCAFFOLD ONLY.
# Ti/Si/C positions are fixed by the matching input, so their internal TSCBOP
# energy is constant and their internal forces do not affect mobile Ag. This
# include omits those fixed-fixed calculations while preserving all Ag forces.
# It is not a flexible-ceramic or physical Ag/MAX interface potential.
pair_style hybrid eam zero 2.40407837 lj/cut 2.40407837
pair_modify pair lj/cut shift yes
pair_coeff * * zero
pair_coeff 4 4 eam potentials/Ag_u3.eam
pair_coeff 1 4 lj/cut 0.5 2.08210927 2.33708864
pair_coeff 2 4 lj/cut 0.5 2.14179034 2.40407837
pair_coeff 3 4 lj/cut 0.5 1.91771227 2.15255924
