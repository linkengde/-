"""Safe dataset checks by default; blind archives only with explicit A runtime flag."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from ase.io import read
PROJECT = Path('research/mace-v12-transfer/periodic_interface_v4')
DRAFT = Path('research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v16_interface_energy')

def require(ok, message):
    if not ok: raise ValueError(message)

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()

def load(path): return json.loads(path.read_text())

def child(root, relative):
    p=(root/relative).resolve(); require(p.is_relative_to(root.resolve()), 'Path escapes root'); return p

def labels(a, energy='REF_energy', forces='REF_forces'):
    require(energy in a.info and forces in a.arrays,'Missing reference fields')
    e=float(a.info[energy]); f=np.asarray(a.arrays[forces],dtype=float)
    require(np.isfinite(e) and f.shape==(len(a),3) and np.isfinite(f).all(),'Nonfinite/invalid reference')
    require(np.isfinite(a.positions).all() and np.isfinite(a.cell.array).all(),'Nonfinite geometry')
    return e,f

def geometry(a,b):
    require(np.array_equal(a.numbers,b.numbers) and np.array_equal(a.pbc,b.pbc),'Wrong elements/order/PBC')
    for k in ['lammps_id','central_pair']:
        require(k in a.arrays and k in b.arrays and np.array_equal(a.arrays[k],b.arrays[k]),'Wrong IDs/marked pair')
    require(len(set(a.arrays['lammps_id']))==len(a),'Duplicate IDs')
    require(np.allclose(a.positions,b.positions,atol=1e-7,rtol=0) and np.allclose(a.cell,b.cell,atol=1e-8,rtol=0),'Wrong geometry/cell')

def checks(repo,archives=False,include_references=False):
    repo=repo.resolve(); project=repo/PROJECT; data=project/'mace_periodic_v16_interface_energy/data'
    pins=load(repo/DRAFT/'input_pins.json')
    for relative,expected in pins['sha256'].items():require(sha(child(repo,relative))==expected,'Pinned input changed: '+relative)
    m=load(data/'dataset_manifest.json'); require(m['split_sizes']=={'train':35,'valid':2,'test':3},'Wrong split counts')
    splits={}
    for s,n in m['split_sizes'].items():
        p=data/(s+'.extxyz');require(sha(p)==m['output_sha256'][s],'Split hash mismatch')
        splits[s]=read(p,index=':');require(len(splits[s])==n,'Split count mismatch')
        for a in splits[s]:labels(a)
    for relative,expected in m['V16_sources_sha256'].items():require(sha(child(repo,relative))==expected,'V16 acquisition source changed: '+relative)
    matrix=np.array([[a.get_chemical_symbols().count(e) for e in ['Ag','C','Si','Ti']] for a in splits['train']])
    require(np.linalg.matrix_rank(matrix)==4 and np.array_equal(matrix,m['energy_composition_matrix']['matrix']),'Composition matrix/rank mismatch')
    require(m['explicit_unresolved_exclusion'] is True,'Approved exclusion missing')
    for s in ['valid','test']:require(sha(data/(s+'.extxyz'))==sha(project/'mace_periodic_v14_interface_energy/data'/(s+'.extxyz')),'Historical split changed')
    # Reuse reviewed geometry screen; no builder execution or data writes.
    import importlib.util
    bp=repo/'research/mace-v12-transfer/coordination/reports/window-b/v15_provenance_aware_builder_revision/build_v15_draft.py'
    require(sha(bp)==m['reviewed_builder_sha256'],'Reviewed builder changed')
    spec=importlib.util.spec_from_file_location('geometry_checks',bp); mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    blind=project/'pbe_interface_v16_blind_holdouts'; manifest=load(blind/'input_manifest.json')
    require({r['label'] for r in manifest['records']}=={f'{x}_registry_holdout_v16_01' for x in ['AgC','AgSi','AgTi']} and len(manifest['records'])==3,'Wrong blind inputs')
    frames=[]
    for r in manifest['records']:
        p=child(blind,r['input']);require(p.parent==blind/'inputs' and sha(p)==r['input_sha256'],'Blind input/hash mismatch')
        a=read(p);require(a.calc is None and not any(k in a.info for k in ['energy','REF_energy','PW_PBE_energy_eV']) and not any(k in a.arrays for k in ['forces','REF_forces','PW_PBE_forces']),'Blind input has labels')
        for s in splits:require(not mod.overlaps(splits[s],[a]),'Blind input leakage: '+s)
        if archives:
            folder=blind/'calculations'/r['label'];v=load(folder/'verification.json');summary=load(folder/'summary.json');progress=load(folder/'progress.json')
            require(v.get('status')=='PASS' and v.get('role')=='v16_blind_local_registry_holdout' and v.get('checks') and all(x is True for x in v['checks'].values()),'Unverified blind archive')
            require(v.get('label')==summary.get('label')==r['label'],'Archive identity mismatch')
            require(v.get('input_sha256')==summary.get('source_sha256')==progress.get('source_sha256')==r['input_sha256'],'Archive input linkage mismatch')
            require(summary.get('scf_converged') is True and progress.get('scf_converged') is True and progress.get('status')=='complete','Incomplete blind SCF')
            for name,h in v['sha256'].items():
                require(Path(name).name==name,'Unsafe archive member');require(sha(folder/name)==h,'Archive hash mismatch: '+name)
            p=folder/(r['label']+'_PW_PBE.extxyz');require(p.name in v['sha256'] and 'summary.json' in v['sha256'] and 'progress.json' in v['sha256'] and 'gpaw.log' in v['sha256'],'Incomplete archive hashes')
            require(summary.get('mpi_ranks')==progress.get('mpi_ranks')==4,'Wrong blind MPI rank count')
            if include_references:
                result=read(p);geometry(a,result);e,f=labels(result,'PW_PBE_energy_eV','PW_PBE_forces')
                require(np.isclose(e,summary['energy_eV_cell'],atol=1e-8,rtol=0),'Archive energy mismatch')
                frames.append((r['label'],result,e,f,sha(p)))
    return {'counts':m['split_sizes'],'rank':4,'blind_archives_checked':archives},splits,frames

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--check-blind-archives',action='store_true');args=ap.parse_args()
    if args.check_blind_archives:
        import subprocess,re
        processes=subprocess.check_output(['ps','-eo','args'],text=True)
        require(not any(re.search(r'(^|/)(mpirun|mpiexec|prterun)(\s|$)',line.strip()) for line in processes.splitlines()),'Active MPI/DFT process; wait for queue completion')
    print(json.dumps(checks(args.repo_root,args.check_blind_archives)[0],indent=2))
