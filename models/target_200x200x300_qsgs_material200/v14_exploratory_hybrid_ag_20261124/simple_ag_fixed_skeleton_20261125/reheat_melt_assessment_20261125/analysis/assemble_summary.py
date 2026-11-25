#!/usr/bin/env python3
"""Merge LAMMPS thermo, explicit-cutoff Q6 rerun, and OVITO PTM results."""
import csv, re, sys
from pathlib import Path
ptm_csv,q6_screen,thermo_screen,out=sys.argv[1:]
ptm={}
with open(ptm_csv,newline='') as f:
    for r in csv.DictReader(f): ptm[int(r['step'])]=r
q6={}
for line in Path(q6_screen).read_text(errors='replace').splitlines():
    x=line.split()
    if len(x)==4 and x[0].isdigit():
        try: q6[int(x[0])]={'Q6_mean_rcut3p8':float(x[1]),'Q6_fraction_gt0p30_rcut3p8':float(x[2])}
        except ValueError: pass
thermo={}
for line in Path(thermo_screen).read_text(errors='replace').splitlines():
    x=line.split()
    if len(x)>=12 and x[0].isdigit():
        try:
            s=int(x[0]); thermo[s]={'time_ps':float(x[1]),'T_Ag_K':float(x[2]),'PE_eV':float(x[3]),'KE_eV':float(x[4]),'Etot_eV':float(x[5]),'pressure_bar':float(x[6]),'Fmax_Ag_eV_per_A':float(x[7]),'Q6_mean_default_cutoff':float(x[8]),'Q6_fraction_gt0p30_default_cutoff':float(x[9]),'MSD_Ag_A2':float(x[10]),'N_total':int(float(x[11]))}
        except ValueError: pass
steps=sorted(set(ptm)&set(q6)&set(thermo))
fields=['step','time_ps','T_Ag_K','PE_eV','KE_eV','Etot_eV','Fmax_Ag_eV_per_A','Q6_mean_rcut3p8','Q6_fraction_gt0p30_rcut3p8','MSD_Ag_A2','N_total','n_Ag','PTM_FCC','PTM_HCP','PTM_BCC','PTM_other','PTM_FCC_HCP_fraction','PTM_other_fraction']
with open(out,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
    for s in steps:
        row={'step':s,**thermo[s],**q6[s]}
        # retain only the explicit-cutoff Q6 values in the final merged table
        row['n_Ag']=int(ptm[s]['n_Ag'])
        for k in ['PTM_FCC','PTM_HCP','PTM_BCC','PTM_other','PTM_FCC_HCP_fraction','PTM_other_fraction']:
            row[k]=ptm[s][k]
        w.writerow({k:row.get(k,'') for k in fields})
print(f'Wrote {len(steps)} joined frames to {out}')
