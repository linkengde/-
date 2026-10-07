import json,hashlib
from pathlib import Path
import numpy as np
from ase.io import read
root=Path(__file__).resolve().parent;rec=json.loads((root/'input_manifest.json').read_text())['records'][0];folder=root/'calculations'/rec['label'];src=root/rec['input'];a=read(src);out=folder/(rec['label']+'_PW_PBE.extxyz');b=read(out);s=json.loads((folder/'summary.json').read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
repo=root.parents[3]; dataset=repo/rec['source_dataset']; assert sha(dataset)==rec['source_dataset_sha256']; original=read(dataset,index=rec['source_frame']);assert np.array_equal(original.positions,a.positions) and np.array_equal(original.cell,a.cell) and np.array_equal(original.numbers,a.numbers) and np.array_equal(original.pbc,a.pbc)
assert s['method']['xc']=='PBE' and s['method']['cutoff_eV']==500 and s['method']['kpts']==[1,1,1] and s['method']['smearing_eV']==.1
assert 'Converged in' in (folder/'gpaw.log').read_text() and np.isfinite(s['free_energy_eV_cell'])
checks={'input_hash':sha(src)==rec['input_sha256'],'positions':np.array_equal(a.positions,b.positions),'cell':np.array_equal(a.cell,b.cell),'pbc':np.array_equal(a.pbc,b.pbc),'symbols':a.get_chemical_symbols()==b.get_chemical_symbols(),'MPI4':s['mpi_ranks']==4,'converged':s['scf_converged'] is True,'energy':b.info['PW_PBE_energy_eV']==s['energy_eV_cell'],'finite':bool(np.isfinite(b.info['PW_PBE_energy_eV']) and np.isfinite(b.arrays['PW_PBE_forces']).all() and b.arrays['PW_PBE_forces'].shape==(len(a),3))}
if not all(checks.values()):raise SystemExit(str(checks))
hashes={p.name:sha(p) for p in folder.iterdir() if p.is_file() and p.suffix!='.gpw' and p.name not in ['verification.json','SHA256SUMS.txt']}
(folder/'verification.json').write_text(json.dumps({'status':'PASS','checks':checks,'sha256':hashes,'source_record':rec},indent=2)+'\n');hashes['verification.json']=sha(folder/'verification.json');(folder/'SHA256SUMS.txt').write_text(''.join(f'{h}  {p}\n' for p,h in hashes.items()));print('archive PASS')
