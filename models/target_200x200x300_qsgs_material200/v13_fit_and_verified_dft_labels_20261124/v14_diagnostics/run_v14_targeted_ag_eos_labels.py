#!/usr/bin/env python3
"""Create a small consistent DFT Ag FCC EOS set for targeted V14 fitting.

The 4x4x4 mesh in the conventional 4-atom cell matches the reciprocal-space
sampling density of the existing 32-atom Ag labels calculated at 2x2x2.
The DFT settings remain exploratory (PW350, sigma=0.10 eV), not universal.
"""
import hashlib,json,time
from pathlib import Path
import numpy as np
from ase.build import bulk
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write
from gpaw import GPAW,FermiDirac,Mixer,PW

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'v14_stage2_p0_ag_dft_diagnostics'/'ag_fcc_eos_k4_equiv'
POINTS=[(4.09,'training'),(3.90,'training'),(3.60,'training'),(3.75,'validation')]
ECUT=350; KPTS=(4,4,4); SIGMA=0.10
TOL={'density':1e-5,'eigenstates':1e-5,'energy':1e-5}

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    OUT.mkdir(parents=True,exist_ok=True); rows=[]
    for a,role in POINTS:
        tag=f'a_{a:.2f}A'; xyz=OUT/f'{tag}.extxyz'; log=OUT/f'{tag}.gpaw.txt'; meta=OUT/f'{tag}.json'
        if meta.exists():
            old=json.loads(meta.read_text())
            if old.get('status')=='complete' and xyz.exists():
                rows.append(old); continue
        atoms=bulk('Ag','fcc',a=a,cubic=True)
        atoms.calc=GPAW(mode=PW(ECUT),xc='PBE',kpts=KPTS,symmetry='off',
          occupations=FermiDirac(SIGMA),mixer=Mixer(beta=0.05,nmaxold=8,weight=100),
          convergence=TOL,maxiter=180,txt=str(log))
        t0=time.time(); row={'candidate_id':f'V14_AG_FCC_EOS_A{int(round(a*100)):04d}',
          'lattice_A':a,'role':role,'natoms':len(atoms),'formula':'Ag4',
          'volume_per_atom_A3':atoms.get_volume()/len(atoms),'nearest_neighbor_A':a/np.sqrt(2),
          'xc':'PBE','code':'GPAW','mode':'plane waves','ecut_eV':ECUT,'kpts':list(KPTS),
          'symmetry':'off','sigma_eV':SIGMA,'training_energy_quantity':'force-consistent finite-smearing free energy',
          'scf_tolerances':TOL,'status':'running'}
        try:
            e_ex=float(atoms.get_potential_energy(force_consistent=False))
            e_free=float(atoms.get_potential_energy(force_consistent=True))
            f=np.asarray(atoms.get_forces(),float); stress=np.asarray(atoms.get_stress(voigt=False),float)
            logtxt=log.read_text(errors='replace') if log.exists() else ''
            conv=bool(getattr(atoms.calc.scf,'converged',False) or ('Converged in ' in logtxt and 'Did not converge!' not in logtxt))
            finite=bool(np.isfinite([e_ex,e_free]).all() and np.isfinite(f).all() and np.isfinite(stress).all())
            row.update({'energy_extrapolated_eV':e_ex,'energy_force_consistent_free_eV':e_free,
              'energy_fit_target_eV':e_free,'forces_eV_A':f.tolist(),'stress_eV_A3':stress.tolist(),
              'force_rmse_vector_norm_eV_A':float(np.sqrt(np.mean(np.sum(f*f,axis=1)))),
              'max_force_norm_eV_A':float(np.linalg.norm(f,axis=1).max()),
              'scf_converged':conv,'finite_energy_force_stress':finite,
              'scf_iterations':int(atoms.calc.get_number_of_iterations()),
              'status':'complete' if conv and finite else 'not_converged'})
            if conv and finite:
                atoms.info['candidate_id']=row['candidate_id']; atoms.info['energy_force_consistent_free_eV']=e_free
                atoms.info['energy_extrapolated_eV']=e_ex; atoms.info['v14_role']=role
                atoms.calc=SinglePointCalculator(atoms,energy=e_free,forces=f,stress=stress.reshape(-1))
                write(xyz,atoms,format='extxyz'); row['extxyz_sha256']=sha(xyz)
        except Exception as exc:
            row.update({'status':'failed','error':f'{type(exc).__name__}: {exc}'})
        row['elapsed_seconds']=round(time.time()-t0,3); row['log_sha256']=sha(log) if log.exists() else None
        meta.write_text(json.dumps(row,indent=2)+'\n'); rows.append(row)
        print(json.dumps({k:row.get(k) for k in ['candidate_id','status','scf_iterations','energy_fit_target_eV','elapsed_seconds']}),flush=True)
    manifest={'protocol':'GPAW/PBE PW350, kpts 4x4x4, sigma 0.10 eV, symmetry off; k-point density matched to V13 Ag32 2x2x2 labels',
      'protocol_status':'exploratory consistent-label set; not claim of universal DFT convergence',
      'rows':[{k:r.get(k) for k in ['candidate_id','lattice_A','role','status','scf_converged','scf_iterations','energy_fit_target_eV','extxyz_sha256','log_sha256']} for r in rows]}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'manifest':str(OUT/'manifest.json'),'complete':sum(r.get('status')=='complete' for r in rows),'total':len(rows)},indent=2))

if __name__=='__main__': main()
