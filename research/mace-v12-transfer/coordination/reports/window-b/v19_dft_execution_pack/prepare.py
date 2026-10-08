import json,hashlib,shutil,ast,subprocess
from pathlib import Path
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parents[6];OUT=Path(__file__).resolve().parent;WB=OUT.parent;BASE=ROOT/'research/mace-v12-transfer/periodic_interface_v4';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
p=WB/'v19_residual_direction_geometry_design/manifest.json';m=json.loads(p.read_text());review_path=BASE/'mace_periodic_v18_clean_core/review_v19_candidates.json';review=json.loads(review_path.read_text());assert review['candidate_manifest_sha256']==sha(p);reviewed={r['candidate_id']:r for r in review['rows']};records=[];(OUT/'inputs').mkdir(exist_ok=True)
for cid in m['recommended']:
 r=next(x for x in m['records'] if x['candidate_id']==cid);src=ROOT/r['path'];assert sha(src)==r['sha256']==reviewed[cid]['sha256'];a=read(src);assert a.calc is None and not any('energy' in k.lower() or 'force' in k.lower() for k in a.info.keys()|a.arrays.keys());assert 'lammps_id' in a.arrays and 'central_pair' in a.arrays
 dest=OUT/'inputs'/src.name;dest.write_bytes(src.read_bytes());b=read(dest);assert np.array_equal(a.positions,b.positions) and np.array_equal(a.cell,b.cell) and np.array_equal(a.arrays['lammps_id'],b.arrays['lammps_id'])
 records.append({'label':cid.replace('_proposal','_acq'),'candidate_id':cid,'input':str(dest.relative_to(OUT)),'input_sha256':sha(dest),'source_path':str(src.relative_to(ROOT)),'role':'training_acquisition_reused_development_family','owner_task':'window-a' if r['sign']>0 else 'window-b','design':r['moved_atoms'],'parent_geometry_sha256':r['parent_geometry_sha256'],'state':'BLOCKED_PENDING_A_INTEGRATION_AND_OWNER_REGISTRATION'})
assert len(records)==6
(OUT/'input_manifest.json').write_text(json.dumps({'launch_enabled':False,'source_manifest_sha256':sha(p),'A_geometry_review_sha256':sha(review_path),'settings':{'xc':'PBE','cutoff_eV':500,'kpts':[1,1,1],'fermi_smearing_eV':.1,'MPI_ranks':4,'BLAS_threads_per_rank':1,'energy':'native extrapolated plus separate free energy'},'records':records},indent=2)+'\n')
# The audited generic frame0 calculator already has no pair/type dependency.
s=(BASE/'pbe_interface_v17_rank_recovery_parallel/run_pw_reference.py').read_text();s='import sys\nraise SystemExit("LAUNCH DISABLED: A must review/integrate and register exact owner jobs; do not enable this report draft in place")\n'+s
s=s.replace('pos:R:3:PW_PBE_forces:R:3','pos:R:3:PW_PBE_forces:R:3:lammps_id:I:1:central_pair:I:1').replace("format(v, '.17g') for v in list(xyz) + list(force)) + '\\n'", "format(v, '.17g') for v in list(xyz) + list(force)) + f' {int(atom_id)} {int(marker)}\\n'")
s=s.replace('for symbol, xyz, force in zip(atoms.get_chemical_symbols(), atoms.positions, forces)','for symbol, xyz, force, atom_id, marker in zip(atoms.get_chemical_symbols(), atoms.positions, forces, atoms.arrays["lammps_id"], atoms.arrays["central_pair"])')
(OUT/'run_pw_reference_draft.py').write_text(s)
for file in ['run_pw_reference_draft.py']:ast.parse((OUT/file).read_text())
print('six exact input hashes, ID/order/cell and label absence PASS; calculator draft syntax PASS')
