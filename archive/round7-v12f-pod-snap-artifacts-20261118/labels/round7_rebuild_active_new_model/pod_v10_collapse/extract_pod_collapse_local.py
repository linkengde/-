#!/workspace/venvs/dft-labels/bin/python
"""Build separate local DFT crops around POD-v10 FIRE collapse contacts."""
from pathlib import Path
import importlib.util
import json
import numpy as np
from ase import Atoms
from ase.io import write

HERE = Path(__file__).resolve().parent
DUMP = Path('/workspace/potential-fit/round7_pod_snap/models/final_v10/pod_minimize_347k_mpi4_restart25/checkpoint.dump')
AUDIT = Path('/workspace/potential-fit/round7_pod_snap/models/final_v10/pod_minimize_347k_conservative/p0_structure_audit.py')
RADIUS = 5.5
CROP_CELL = 22.0
CENTERS = [
    (76905, 76907, 'Ti-Si collapse center on Ti'),
    (76907, 76905, 'Ti-Si collapse center on Si'),
    (140, 2748, 'C-C collapse center on C 140'),
    (2748, 140, 'C-C collapse center on C 2748'),
]
SYMBOL = {1:'Ag',2:'C',3:'Si',4:'Ti'}

spec=importlib.util.spec_from_file_location('p0audit',AUDIT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
step, ids, typ, xyz, bounds = mod.read_frame(DUMP,-1)
if step != 50:
    raise RuntimeError(f'Expected last saved frame at step 50, got {step}')
id_to_i={int(atom_id):i for i,atom_id in enumerate(ids)}
lx=bounds[0,1]-bounds[0,0]
ly=bounds[1,1]-bounds[1,0]
frames=[]
records=[]
for center_id, partner_id, reason in CENTERS:
    ci=id_to_i[center_id]
    pi=id_to_i[partner_id]
    delta=xyz-xyz[ci]
    delta[:,0]-=lx*np.rint(delta[:,0]/lx)
    delta[:,1]-=ly*np.rint(delta[:,1]/ly)
    dist=np.linalg.norm(delta,axis=1)
    keep=np.flatnonzero(dist<=RADIUS)
    pos=delta[keep]+CROP_CELL/2
    local=Atoms([SYMBOL[int(typ[i])] for i in keep],positions=pos,
                cell=np.eye(3)*CROP_CELL,pbc=True)
    local.new_array('source_atom_id',ids[keep].astype(np.int32))
    local.info.update({
        'config_type':f'POD_v10_FIRE_step50_{reason.replace(" ","_")}',
        'geometry_source':str(DUMP),'source_step':int(step),
        'center_atom_id':int(center_id),'center_species':SYMBOL[int(typ[ci])],
        'partner_atom_id':int(partner_id),'partner_species':SYMBOL[int(typ[pi])],
        'center_partner_distance_A':float(dist[pi]),'crop_radius_A':RADIUS,
        'crop_cell_A':CROP_CELL,'source_pbc':'x/y periodic, z open; local crop embedded in vacuum box',
        'label_role':'training candidate only; not an independent blind-test structure',
        'boundary_note':'Finite local crop; use center and contact-partner forces for exploratory fitting, downweight edge atoms; raw energy cannot be directly mixed with full-periodic training energies.',
    })
    frames.append(local)
    records.append({'center_atom_id':center_id,'center_species':SYMBOL[int(typ[ci])],
                    'partner_atom_id':partner_id,'partner_species':SYMBOL[int(typ[pi])],
                    'pair_distance_A':float(dist[pi]),'crop_atoms':len(local),
                    'formula':local.get_chemical_formula(),'reason':reason})
out=HERE/'candidate_geometries.extxyz'
write(out,frames,format='extxyz')
(HERE/'extraction_summary.json').write_text(json.dumps({'source_dump':str(DUMP),'source_step':step,
    'radius_A':RADIUS,'cell_A':CROP_CELL,'role':'training-candidate; separate from blind test',
    'frames':records},indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'output':str(out),'frames':records},indent=2,ensure_ascii=False))
