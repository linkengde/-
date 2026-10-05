#!/usr/bin/env python3
"""A review-only compact archive verifier; no DFT import, execution or job claim."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
from ase.io import read

def require(ok,message):
    if not ok:raise ValueError(message)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--run-dir',required=True,type=Path);ap.add_argument('--archive-dir',required=True,type=Path);ap.add_argument('--label',required=True);args=ap.parse_args()
    repo=args.repo_root.resolve();pack=repo/'research/mace-v12-transfer/coordination/reports/window-b/next_acquisition_execution_pack';run=args.run_dir.resolve();dest=args.archive_dir.resolve()
    require(dest.is_relative_to(pack) and dest!=pack,'Archive must be isolated report child')
    require(not dest.exists(),'Archive exists: refuse overwrite even empty directories')
    records=json.loads((pack/'input_manifest.json').read_text())['records'];record=next((x for x in records if x['label']==args.label),None);require(record is not None,'Unknown planned label')
    inp=(pack/record['input']).resolve();require(inp.parent==pack/'inputs' and sha(inp)==record['input_sha256'],'Input hash/path changed')
    require(not (run/'active_job.json').exists(),'Active/failed run marker exists; preserve results and investigate')
    summary=json.loads((run/'summary.json').read_text());progress=json.loads((run/'progress.json').read_text());log=(run/'gpaw.log').read_text()
    require(summary['label']==progress['label']==record['label'] and summary['source_sha256']==progress['source_sha256']==record['input_sha256'],'Output linkage mismatch')
    require(summary['scf_converged'] is True and progress.get('scf_converged') is True and progress['status']=='complete','Not converged/complete')
    require(summary['mpi_ranks']==progress['mpi_ranks']==4 and f"Converged in {summary['scf_iterations']} steps" in log,'MPI/convergence log mismatch')
    method=summary['method'];require(method['xc']=='PBE' and method['basis']=='plane wave' and method['cutoff_eV']==500 and method['kpts']==[1,1,1] and method['smearing_eV']==.1,'Unexpected method')
    source=read(inp);path=run/(args.label+'_PW_PBE.extxyz');a=read(path)
    require(np.array_equal(a.numbers,source.numbers) and np.array_equal(a.pbc,source.pbc),'Element/order/PBC mismatch')
    for key in ['lammps_id','central_pair']:require(np.array_equal(a.arrays[key],source.arrays[key]),'ID/pair mismatch')
    require(len(set(a.arrays['lammps_id']))==len(a),'Duplicate IDs')
    require(np.allclose(a.positions,source.positions,atol=1e-12,rtol=0) and np.allclose(a.cell,source.cell,atol=1e-12,rtol=0),'Geometry/cell mismatch')
    energy=float(a.info['PW_PBE_energy_eV']);forces=a.arrays['PW_PBE_forces'];require(np.isfinite(energy) and forces.shape==(len(a),3) and np.isfinite(forces).all(),'Invalid labels')
    require(np.isclose(energy,summary['energy_eV_cell'],atol=1e-10,rtol=0),'Chosen energy mismatch')
    free=summary.get('free_energy_eV_cell');require(free is None or (np.isfinite(free) and 'PW_PBE_free_energy_eV' in a.info and np.isclose(free,a.info['PW_PBE_free_energy_eV'],atol=1e-10,rtol=0)),'Free energy mismatch')
    require('force_consistent=False' in summary['energy_convention'],'Unknown chosen convention')
    names=[path.name,'summary.json','progress.json','gpaw.log'];hashes={n:sha(run/n) for n in names}
    dest.parent.mkdir(parents=True,exist_ok=True);staging=dest.parent/(dest.name+'.staging');require(not staging.exists(),'Existing staging: preserve diagnostics');staging.mkdir()
    try:
        for name in names:shutil.copy2(run/name,staging/name);require(sha(staging/name)==hashes[name],'Copy hash mismatch')
        verify={'status':'PASS','label':args.label,'role':record['role'],'input_sha256':record['input_sha256'],'sha256':hashes,'geometry_IDs_finite_energy_forces_convergence_checks':True,'free_energy_recorded':free is not None,'REF_convention':'PW_PBE_energy_eV only; do not substitute free energy','no_gpw_archived':True}
        (staging/'verification.json').write_text(json.dumps(verify,indent=2)+'\n');staging.rename(dest)
    except Exception as e:
        (staging/'diagnostics_failure.json').write_text(json.dumps({'error':str(e)})+'\n');raise
if __name__=='__main__':main()
