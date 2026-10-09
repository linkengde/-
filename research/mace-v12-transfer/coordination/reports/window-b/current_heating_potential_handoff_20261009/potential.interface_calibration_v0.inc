# CALIBRATION VERSION 0 - NOT YET A PRODUCTION POTENTIAL
# Atom types: 1 Ti, 2 Si, 3 C, 4 Ag

pair_style hybrid eam tersoff morse 8.0

# Phase-internal interactions
pair_coeff 4 4 eam potentials/Ag_u3.eam
pair_coeff * * tersoff potentials/TSCBOP.tersoff Ti Si C NULL

# Cross-interface Morse starting values fitted to a single-plane HCP surrogate
# of the reported Ag(111)/Ti3SiC2(0001) DFT adhesion targets.
# pair_coeff i j morse D0[eV] alpha[1/A] r0[A] cutoff[A]
pair_coeff 1 4 morse 0.94814194 2.00000000 2.33708864 8.0
pair_coeff 2 4 morse 0.93944948 2.00000000 2.40407837 8.0
pair_coeff 3 4 morse 1.35077230 2.00000000 2.15255924 8.0

