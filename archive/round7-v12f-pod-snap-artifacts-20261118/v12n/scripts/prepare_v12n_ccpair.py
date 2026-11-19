#!/workspace/venvs/dft-labels/bin/python
"""Create isolated v12n POD/SNAP candidate from v12m plus two DFT C-C labels."""
from __future__ import annotations
import hashlib,json,shutil,subprocess
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read,write
HERE=Path(__file__).resolve().parent
V12M=HERE/'v12m_agctriad_candidate'; OUT=HERE/'v12n_ccpair_candidate'
ACTIVE=HERE.parent/'target_200x200x300_qsgs_material200'/'active_tic_hotspot_expansion'
LABELS=[('cc_pair_v12m_t26fs_step16_pw350_gamma',[128985,134363]),('cc_pair_v12m_t28fs_step18_pw350_gamma',[128985,134363])]
V10=HERE.parent/'final_v10'/'POD'; LMP=Path('/workspace/software/lammps-pod-snap-clean/build/lmp')
ELEMENTS=['Ag','C','Si','Ti']; MASS={'Ag':107.8682,'C':12.011,'Si':28.085,'Ti':47.867}
def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def run_v10_pod(atoms,work):
    work.mkdir(parents=True,exist_ok=True); data=work/'local.data'
    write(data,atoms,format='lammps-data',atom_style='atomic',units='metal',specorder=ELEMENTS)
    masses='\n'.join(f'mass {i} {MASS[s]}' for i,s in enumerate(ELEMENTS,1))
    inp=work/'forces.in'; inp.write_text(f'''clear\nunits metal\natom_style atomic\nboundary p p p\nread_data {data}\n{masses}\npair_style pod\npair_coeff * * {V10/'AgTiSiC_POD.pod'} {V10/'pod_coefficients.pod'} Ag C Si Ti\nneighbor 2.0 bin\nneigh_modify delay 0 every 1 check yes\nrun 0\nwrite_dump all custom forces.dump id type x y z fx fy fz modify sort id\n''')
    subprocess.run([str(LMP),'-in',str(inp),'-log',str(work/'lammps.log'),'-screen','none'],cwd=work,check=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    lines=(work/'forces.dump').read_text().splitlines(); h=next(i for i,x in enumerate(lines) if x.startswith('ITEM: ATOMS ')); cols=lines[h].split()[2:]
    a=np.loadtxt(lines[h+1:],usecols=[cols.index(k) for k in ('id','fx','fy','fz')],dtype=float,ndmin=2); a=a[np.argsort(a[:,0])]
    if len(a)!=len(atoms) or not np.isfinite(a[:,1:]).all():raise RuntimeError('Invalid v10 anchor forces')
    return a[:,1:4]
def main():
    if OUT.exists():raise FileExistsError(OUT)
    for p in (V12M/'dataset/local_train/local_dft_targets.xyz',V12M/'POD/AgTiSiC_POD_v12m.pod',V12M/'POD/fitpod.data',V12M/'SNAP/fitsnap.in',V10/'AgTiSiC_POD.pod',V10/'pod_coefficients.pod',LMP):
        if not p.exists():raise FileNotFoundError(p)
    records=[]
    for stem,target_ids in LABELS:
        xyz=ACTIVE/'labels'/f'{stem}.extxyz'; js=xyz.with_suffix('.json'); log=xyz.with_suffix('.txt')
        if not all(p.exists() for p in (xyz,js,log)):raise FileNotFoundError(xyz)
        rec=json.loads(js.read_text())
        if rec.get('status')!='CONVERGED' or 'Converged in ' not in log.read_text(errors='replace'):raise RuntimeError(f'{stem}: DFT label is not SCF-converged')
        if rec.get('center_atom_ids')!=target_ids or rec.get('DFT_target_rows')!=2:raise RuntimeError(f'{stem}: C-C target IDs/rows mismatch')
        records.append((stem,target_ids,xyz,js,log,rec))
    shutil.copytree(V12M/'dataset/train',OUT/'dataset/train',ignore=shutil.ignore_patterns('local_dft_targets.xyz'))
    shutil.copytree(V12M/'dataset/validation',OUT/'dataset/validation')
    shutil.copytree(V12M/'SNAP/data',OUT/'SNAP/data',ignore=shutil.ignore_patterns('Local.xyz'))
    (OUT/'dataset/local_train').mkdir(parents=True); (OUT/'POD').mkdir(); (OUT/'SNAP').mkdir(exist_ok=True)
    for src,dst in [(V12M/'holdout_si_hotspot_104359.extxyz',OUT/'holdout_si_hotspot_104359.extxyz'),(V12M/'holdout_label_evaluation/POD_v12f_si01/local.data',OUT/'holdout_label_evaluation/POD_v12f_si01/local.data')]:
        dst.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(src,dst)
    shutil.copy2(V12M/'POD/AgTiSiC_POD_v12m.pod',OUT/'POD/AgTiSiC_POD_v12n.pod'); shutil.copy2(V12M/'SNAP/fitsnap.in',OUT/'SNAP/fitsnap.in')
    local=read(V12M/'dataset/local_train/local_dft_targets.xyz',index=':')
    if len(local)!=22:raise RuntimeError(f'Expected 22 inherited local frames, got {len(local)}')
    added=[]
    for stem,target_ids,xyz,js,log,rec in records:
        dft=read(xyz); ids=np.asarray(dft.arrays['source_atom_id'],int); target=np.asarray(dft.arrays['is_DFT_force_target'],bool); ref=np.asarray(dft.arrays['REF_forces'],float)
        if target.sum()!=2 or sorted(ids[target].tolist())!=sorted(target_ids) or not np.isfinite(ref).all():raise RuntimeError(f'{stem}: bad DFT target mapping/data')
        base=run_v10_pod(dft,OUT/'local_pod_v10_forces'/stem); fit=base.copy(); fit[target]=ref[target]
        rich=Atoms(symbols=dft.get_chemical_symbols(),positions=dft.positions,cell=dft.cell,pbc=dft.pbc)
        rich.new_array('forces',fit); rich.new_array('source_atom_id',ids.astype(np.int32)); rich.new_array('DFT_forces',ref); rich.new_array('POD_v10_forces',base); rich.new_array('is_DFT_force_target',target.astype(np.int8))
        rich.info.update({'energy':float(dft.info['REF_energy']),'config_type':f'POD_v12n_full200_CC_{stem}','center_atom_ids':json.dumps(target_ids),'DFT_method':str(dft.info['DFT_method']),'DFT_label_sha256':sha256(xyz),'role':'training_candidate; Si01 center 104359 remains held out'})
        local.append(rich); shutil.copy2(xyz,OUT/f'{stem}.extxyz'); shutil.copy2(js,OUT/f'{stem}.json'); shutil.copy2(log,OUT/f'{stem}.txt')
        ids_index={int(a):i for i,a in enumerate(ids)}
        added.append({'label':stem,'sha256':sha256(xyz),'target_atom_ids':target_ids,'target_rows':2,'distance_A':rec['center_pair_distance_A'],'force_norms_eV_A':[float(np.linalg.norm(ref[ids_index[a]])) for a in target_ids],'source_model_errors_eV_A':rec['center_POD_v12m_force_errors_eV_A'],'SCF_iterations':rec['iterations']})
    if len(local)!=24:raise RuntimeError(f'Expected 24 local frames total, got {len(local)}')
    write(OUT/'dataset/local_train/local_dft_targets.xyz',local,format='extxyz')
    fitframes=[]
    for src in local:
        forces=np.asarray(src.arrays['forces'] if 'forces' in src.arrays else src.get_forces(apply_constraint=False),float)
        clean=Atoms(symbols=src.get_chemical_symbols(),positions=src.positions,cell=src.cell,pbc=src.pbc); clean.new_array('forces',forces.copy()); clean.info['energy']=float(src.info['energy'] if 'energy' in src.info else src.get_potential_energy()); clean.info['config_type']=str(src.info['config_type']); fitframes.append(clean)
    for path in [OUT/'dataset/train/local_dft_targets.xyz',OUT/'SNAP/data/Local.xyz']:
        write(path,fitframes,format='extxyz'); check=read(path,index=':')
        if len(check)!=24:raise RuntimeError(f'Expected 24 local training frames in {path}')
        for i,(a,b) in enumerate(zip(fitframes,check)):
            if not np.array_equal(a.numbers,b.numbers) or not np.allclose(a.positions,b.positions,atol=0,rtol=0) or not np.allclose(a.arrays['forces'],b.get_forces(),atol=1e-12,rtol=0):raise RuntimeError(f'Round-trip mismatch {path} frame {i}')
    cfg=(V12M/'POD/fitpod.data').read_text().replace(str(V12M),str(OUT)).replace('basename_for_output_files pod_v12m','basename_for_output_files pod_v12n'); (OUT/'POD/fitpod.data').write_text(cfg); (OUT/'POD/fitpod.in').write_text('fitpod AgTiSiC_POD_v12n.pod fitpod.data\n')
    snap=(OUT/'SNAP/fitsnap.in').read_text().replace('SNAP_v12m_metrics.md','SNAP_v12n_metrics.md').replace('AgTiSiC_SNAP_v12m','AgTiSiC_SNAP_v12n'); (OUT/'SNAP/fitsnap.in').write_text(snap)
    manifest=json.loads((V12M/'local_training_manifest.json').read_text()); manifest.update({'status':'prepared; waiting for independent fit and validation','candidate':'v12n_full200_POD_SNAP_C-C_cold-stability_refit','parent_candidate':'v12m_agctriad_candidate','local_force_frames':24,'baseline_train_configurations':128,'baseline_validation_configurations':12,'si01_holdout_center_atom_id':104359,'holdout_excluded_from_training':True,'total_fitpod_training_configurations':152,'new_cc_cold_hotspot_labels':added,'total_DFT_target_force_rows':21,'note':'Preserves 128 train / 12 validation split and Si01 center holdout. Only central C force rows in the two new crops use DFT; crop-edge forces use v10 POD anchors.'})
    manifest['source_labels']=list(manifest.get('source_labels',[]))+[str(x[2]) for x in records]
    (OUT/'local_training_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    # Reuse the established evaluator but bind it to the v12n output names.
    ev=(V12M/'evaluate_v12m_si01.py').read_text().replace('v12m','v12n'); (OUT/'evaluate_v12n_si01.py').write_text(ev)
    fitsh=OUT/'run_v12n_fits.sh'; fitsh.write_text('''#!/usr/bin/env bash\nset -euo pipefail\nroot="$(cd "$(dirname "$0")" && pwd)"\n(cd "$root/POD" && "/workspace/software/lammps-pod-snap-clean/build/lmp" -in fitpod.in -log fitpod_v12n.log)\ntest -s "$root/POD/pod_v12n_coefficients.pod"\ntest -s "$root/POD/pod_v12n_test_errors.pod"\n(cd "$root/SNAP" && "/workspace/venvs/fitsnap/bin/python" /workspace/potential-fit/round7_pod_snap/run_fitsnap_py312.py --overwrite fitsnap.in > fitsnap_v12n.log 2>&1)\ntest -s "$root/SNAP/AgTiSiC_SNAP_v12n.snapcoeff"\ntest -s "$root/SNAP/AgTiSiC_SNAP_v12n.snapparam"\ntest -s "$root/SNAP/SNAP_v12n_metrics.md"\n'''); fitsh.chmod(0o755)
    (OUT/'README_v12n.md').write_text('''# v12n POD/SNAP candidate\n\nParent: v12m. This candidate preserves the existing 128 training / 12 validation split and held-out Si01 center 104359. It adds two GPAW PBE/PW350 Gamma DFT force labels for the same C-C pair at 1.386 and 1.114 Å during a cold-run collapse precursor. Only the central two C force rows per crop are DFT targets; edge rows use v10 POD force anchors.\n\nThis is a training candidate only. It must pass the independent Si01 holdout, full-model static force checks, and full-model cold stability before any heating trial.\n''')
    print(json.dumps({'status':'PREPARED','candidate':str(OUT),'train_base':128,'validation':12,'local_frames':24,'fit_configurations':152,'new_DFT_target_rows':4,'total_DFT_target_force_rows':21,'holdout_center':104359,'labels':added},indent=2))
if __name__=='__main__':main()
