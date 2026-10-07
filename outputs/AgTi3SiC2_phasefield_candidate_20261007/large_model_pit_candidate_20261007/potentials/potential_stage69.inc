# Stage69 collision-ready candidate potential; no physical run is implied.
# Baseline: Ag EAM + Ti/Si/C TSCBOP + calibrated Ag-interface Morse.
# ZBL is overlaid only in the 0.8-1.5 A short-range window.
pair_style hybrid/overlay eam tersoff morse 8.0 zbl 0.8 1.5
pair_coeff 4 4 eam potentials/Ag_u3.eam
pair_coeff * * tersoff potentials/TSCBOP.tersoff Ti Si C NULL
pair_coeff 1 4 morse 0.94814194 2.00000000 2.33708864 8.0
pair_coeff 2 4 morse 0.93944948 2.00000000 2.40407837 8.0
pair_coeff 3 4 morse 1.35077230 2.00000000 2.15255924 8.0
pair_coeff 1 4 zbl 22 47
pair_coeff 2 4 zbl 14 47
pair_coeff 3 4 zbl 6 47
pair_coeff 4 4 zbl 47 47
