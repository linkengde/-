from pathlib import Path
import json
import numpy as np
from ase import Atoms
from ase.io import write

HERE=Path(__file__).resolve().parents[1]
OUT=HERE/'active_tic_hotspot_expansion'
PAIR=(62864,62871)  # Ti-C pair: 2.079 A initially, 1.683 A after 10 fs cold v12g
RADIUS=5.5
CELL=22.0
MODEL_TYPE_TO_ELEMENT={1:'Ag',2:'C',3:'Si',4:'Ti'}
SOURCE_TYPE_TO_ELEMENT={1:'Ti',2:'Si',3:'C',4:'Ag'}

def dump(path):
    lines=Path(path).read_text().splitlines()
    atom_header=next(i for i,s in enumerate(lines) if s.startswith('ITEM: ATOMS'))
    columns=lines[atom_header].split()[2:]
    rows=np.loadtxt(lines[atom_header+1:],dtype=float)
    return columns,rows

def load_data():
    lines=(HERE/'Ag_Ti3SiC2_bicontinuous_200x200x300.data').read_text().splitlines()
    start=next(i for i,s in enumerate(lines) if s.startswith('Atoms'))+1
    rows=[]
    for line in lines[start+1:]:
        if not line.strip():
            continue
        row=line.split()
        if len(row)<5: break
        rows.append([float(x) for x in row[:5]])
        if len(rows)==458248: break
    a=np.asarray(rows)
    return {int(r[0]):(SOURCE_TYPE_TO_ELEMENT[int(r[1])],r[2:5]) for r in a}

source=load_data()
configs=[
 ('initial', HERE/'static_v12g_initial_458k/forces.dump'),
 ('cold_10fs', HERE/'cold_v12g_initial_458k_1fs/end_0p01ps.dump'),
]
summary={'pair_ids':list(PAIR),'source_pair_elements':[source[i][0] for i in PAIR], 'crop_radius_A':RADIUS,'crop_cell_A':CELL,'candidates':[]}
for name,path in configs:
    cols,rows=dump(path)
    ix={k:cols.index(k) for k in ['id','type','x','y','z','fx','fy','fz']}
    data={int(r[ix['id']]):r for r in rows}
    if any(i not in data for i in PAIR): raise RuntimeError(f'{name}: missing pair ids')
    p0=np.array([data[PAIR[0]][[ix['x'],ix['y'],ix['z']]],data[PAIR[1]][[ix['x'],ix['y'],ix['z']]]])
    delta=p0[1]-p0[0]; delta[:2]-=200*np.rint(delta[:2]/200)
    center=p0[0]+0.5*delta
    ids=np.asarray(sorted(data),dtype=np.int64)
    coords=np.asarray([data[int(i)][[ix['x'],ix['y'],ix['z']]] for i in ids])
    d=coords-center; d[:,:2]-=200*np.rint(d[:,:2]/200)
    mask=np.linalg.norm(d,axis=1)<=RADIUS
    kept=ids[mask]; positions=d[mask]+CELL/2
    symbols=[MODEL_TYPE_TO_ELEMENT[int(data[int(i)][ix['type']])] for i in kept]
    forces=np.asarray([data[int(i)][[ix['fx'],ix['fy'],ix['fz']] ] for i in kept])
    target=np.isin(kept,PAIR).astype(np.int8)
    atoms=Atoms(symbols=symbols,positions=positions,cell=[CELL,CELL,CELL],pbc=True)
    atoms.new_array('source_atom_id',kept.astype(np.int32))
    atoms.new_array('POD_v12g_forces',forces)
    atoms.new_array('is_DFT_force_target',target)
    atoms.info.update({'config_type':f'full200_tic_pair_{name}', 'source_timestep':0 if name=='initial' else 10,
        'center_atom_id':int(PAIR[0]),'center_species':source[PAIR[0]][0], 'partner_atom_id':int(PAIR[1]),'partner_species':source[PAIR[1]][0],
        'center_partner_distance_A':float(np.linalg.norm(delta)), 'crop_radius_A':RADIUS,'crop_cell_A':CELL,
        'geometry_source':str(path),'source_pbc':'x/y periodic, z open; crop embedded in cubic vacuum box',
        'label_role':'training candidate; SI01 center 104359 remains held out',
        'DFT_edge_force_note':'Only Ti and C pair center forces are used for training; local crop edge forces excluded.',
        'POD_v12g_force_pair_eV_A':json.dumps(forces[np.where(target==1)[0]].tolist())})
    out=OUT/'candidates'/f'tic_pair_{name}.extxyz'
    write(out,atoms,format='extxyz')
    # Starting and post-10fs relative acceleration projection using full-system model forces.
    idx={int(v):i for i,v in enumerate(kept)}
    i,j=idx[PAIR[0]],idx[PAIR[1]]
    fti=forces[i] if symbols[i]=='Ti' else forces[j]
    fc=forces[i] if symbols[i]=='C' else forces[j]
    # displacement vector from Ti to C in this crop
    rtic=atoms.positions[j]-atoms.positions[i] if symbols[i]=='Ti' else atoms.positions[i]-atoms.positions[j]
    rhat=rtic/np.linalg.norm(rtic)
    relproj=float(np.dot(fc/12.011-fti/47.867,rhat))
    rec={'state':name,'natoms':len(atoms),'pair_distance_A':float(np.linalg.norm(delta)),
         'target_source_ids':list(PAIR),'target_symbols':[source[i][0] for i in PAIR],
         'target_forces_v12g_eV_A':{str(PAIR[0]):forces[idx[PAIR[0]]].tolist(),str(PAIR[1]):forces[idx[PAIR[1]]].tolist()},
         'relative_acceleration_projection_mass_weighted_eV_A_amu_toward_separation':relproj,
         'candidate_extxyz':str(out)}
    summary['candidates'].append(rec)
(OUT/'candidate_manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
