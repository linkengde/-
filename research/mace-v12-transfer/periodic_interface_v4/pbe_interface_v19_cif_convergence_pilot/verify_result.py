"""Read-only compact-result verification draft; no GPAW imports or launch."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from ase.io import read
root=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('label');ap.add_argument('result_dir',type=Path);args=ap.parse_args();r=next(x for x in json.loads((root/'input_manifest.json').read_text())['records'] if x['label']==args.label);folder=args.result_dir;src=root/r['input'];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();a=read(src);p=folder/(args.label+'_PW_PBE.extxyz');b=read(p);s=json.loads((folder/'summary.json').read_text());progress=json.loads((folder/'progress.json').read_text());log=(folder/'gpaw.log').read_text()
checks={'input_hash':sha(src)==r['input_sha256']==s['source_sha256'],'symbol_order':a.get_chemical_symbols()==b.get_chemical_symbols(),'positions':np.array_equal(a.positions,b.positions),'cell':np.array_equal(a.cell,b.cell),'pbc':np.array_equal(a.pbc,b.pbc),'ids':np.array_equal(a.arrays['lammps_id'],b.arrays['lammps_id']),'pair_markers':np.array_equal(a.arrays['central_pair'],b.arrays['central_pair']),'convergence':s['scf_converged'] is True and progress['status']=='complete' and f"Converged in {s['scf_iterations']} steps" in log,'mpi4':s['mpi_ranks']==4,'energy_exact':b.info['PW_PBE_energy_eV']==s['energy_eV_cell'],'finite':np.isfinite(s['energy_eV_cell']) and np.isfinite(s['free_energy_eV_cell']) and b.arrays['PW_PBE_forces'].shape==(len(a),3) and np.isfinite(b.arrays['PW_PBE_forces']).all(),'method':s['method']['xc']=='PBE' and s['method']['cutoff_eV']==500 and s['method']['kpts']==r['kpts'] and s['method']['smearing_eV']==.1}
checks={k:bool(v) for k,v in checks.items()}
if not all(checks.values()):raise SystemExit(json.dumps(checks))
print(json.dumps({'status':'PASS','checks':checks,'owner_task':r['owner_task'],'sha256':{p.name:sha(p) for p in folder.iterdir() if p.is_file() and p.suffix!='.gpw'},'input_sha256':r['input_sha256'],'large_gpw_archived':False},indent=2))
