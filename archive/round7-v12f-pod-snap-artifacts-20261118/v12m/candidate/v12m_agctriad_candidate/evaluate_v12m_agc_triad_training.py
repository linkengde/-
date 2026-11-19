#!/workspace/venvs/dft-labels/bin/python
"""Compare fitted v12m forces to the three DFT targets in its Ag-C-C training crop."""
from __future__ import annotations
import json, subprocess
from pathlib import Path
import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
ACTIVE = Path('/workspace/potential-fit/round7_pod_snap/models/target_200x200x300_qsgs_material200/active_tic_hotspot_expansion')
STEM = 'agc_triad_v12l_t28fs_step18_pw350_gamma'
LABEL = ACTIVE / 'labels' / f'{STEM}.extxyz'
LMP = Path('/workspace/software/lammps-pod-snap-clean/build/lmp')
ELEMENTS = ['Ag', 'C', 'Si', 'Ti']
MASS = 'mass 1 107.8682\nmass 2 12.011\nmass 3 28.085\nmass 4 47.867'
TARGETS = [414900, 198354, 203868]
PAIR = {
    'POD_v12m': f"pair_style pod\npair_coeff * * {ROOT/'POD/AgTiSiC_POD_v12m.pod'} {ROOT/'POD/pod_v12m_coefficients.pod'} Ag C Si Ti",
    'SNAP_v12m': f"pair_style snap\npair_coeff * * {ROOT/'SNAP/AgTiSiC_SNAP_v12m.snapcoeff'} {ROOT/'SNAP/AgTiSiC_SNAP_v12m.snapparam'} Ag C Si Ti",
}
atoms = read(LABEL)
ids = np.asarray(atoms.arrays['source_atom_id'], dtype=int)
target = np.asarray(atoms.arrays['is_DFT_force_target'], dtype=bool)
ref = np.asarray(atoms.arrays['REF_forces'], dtype=float)
if ids[target].tolist() != [198354, 203868, 414900] and sorted(ids[target].tolist()) != sorted(TARGETS):
    raise RuntimeError('DFT target atom IDs do not match the Ag-C-C label')
index = {int(atom_id): i for i, atom_id in enumerate(ids)}
results=[]
for name, pair in PAIR.items():
    work=ROOT/'agc_triad_training_diagnostics'/name
    work.mkdir(parents=True, exist_ok=True)
    data=work/'local.data'
    write(data, atoms, format='lammps-data', atom_style='atomic', units='metal', specorder=ELEMENTS)
    dump=work/'forces.dump'
    inp=work/'run.in'
    inp.write_text(f'''clear\nunits metal\natom_style atomic\nboundary p p p\nread_data {data}\n{MASS}\n{pair}\nneighbor 2.0 bin\nneigh_modify delay 0 every 1 check yes\nrun 0\nwrite_dump all custom {dump} id type fx fy fz modify sort id\n''')
    subprocess.run([str(LMP), '-in', str(inp), '-log', str(work/'log.lammps'), '-screen', 'none'], cwd=work, check=True)
    lines=dump.read_text().splitlines(); h=next(i for i,x in enumerate(lines) if x.startswith('ITEM: ATOMS ')); cols=lines[h].split()[2:]
    ix={k:cols.index(k) for k in ('id','fx','fy','fz')}
    table={}
    for line in lines[h+1:]:
        vals=line.split(); table[int(float(vals[ix['id']]))]=np.array([float(vals[ix[k]]) for k in ('fx','fy','fz')])
    model=np.array([table[index[i]+1] for i in TARGETS])
    dft=np.array([ref[index[i]] for i in TARGETS])
    per=np.linalg.norm(model-dft,axis=1)
    results.append({'model':name,'target_ids':TARGETS,'DFT_force_vectors_eV_A':dft.tolist(),'model_force_vectors_eV_A':model.tolist(),'per_atom_vector_error_eV_A':per.tolist(),'target_mean_vector_error_eV_A':float(per.mean()),'is_training_diagnostic':True})
out={'status':'EVALUATED_TRAINING_DIAGNOSTIC','DFT_label':str(LABEL),'results':results}
(ROOT/'agc_triad_training_diagnostics/summary.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
