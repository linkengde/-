#!/workspace/venvs/gpaw/bin/python
"""Label one V13 Ag-liquid pre-collapse frame for V14 active-learning review.
This is a candidate training label, not an independent validation sample.
"""
import hashlib, json, time, traceback
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
from gpaw import GPAW, FermiDirac, Mixer, PW

BASE=Path(__file__).resolve().parent
OUT=BASE/'v14_stage2_p0_ag_dft_diagnostics'/'ag_liquid_precollapse'
MANIFEST=OUT/'candidate_manifest.json'

m=json.loads(MANIFEST.read_text())
geom=Path(m['geometry_file'])
log=OUT/'gpaw_scf.txt'; result=OUT/'label.json'; label=OUT/'label.extxyz'
a=read(geom,format='extxyz')
rec={**m,'dft_code':'GPAW','xc':'PBE','ecut_eV':350,'kpts':[2,2,2],
     'symmetry':'off','sigma_eV':0.10,
     'energy_fit_target':'force-consistent finite-smearing free energy',
     'scf_tolerances':{'density':1e-5,'eigenstates':1e-5,'energy':1e-5},
     'max_scf_iterations':250,'run_status':'running'}
t0=time.time()
try:
    a.calc=GPAW(mode=PW(350),xc='PBE',kpts=(2,2,2),symmetry='off',
        occupations=FermiDirac(0.10),mixer=Mixer(beta=0.05,nmaxold=8,weight=100),
        convergence={'density':1e-5,'eigenstates':1e-5,'energy':1e-5},
        maxiter=250,txt=str(log))
    e_extrap=float(a.get_potential_energy(force_consistent=False))
    e_free=float(a.get_potential_energy(force_consistent=True))
    f=np.asarray(a.get_forces(),float)
    stress=np.asarray(a.get_stress(voigt=False),float)
    txt=log.read_text(errors='replace').replace('\x00','') if log.exists() else ''
    conv=bool(getattr(getattr(a.calc,'scf',None),'converged',False) or ('Converged in ' in txt and 'Did not converge!' not in txt))
    finite=bool(np.isfinite(e_extrap) and np.isfinite(e_free) and np.isfinite(f).all() and np.isfinite(stress).all())
    rec.update({'energy_extrapolated_eV':e_extrap,'energy_force_consistent_free_eV':e_free,
        'energy_fit_target_eV':e_free,'free_energy_eV_per_atom':e_free/len(a),
        'extrapolated_energy_eV_per_atom':e_extrap/len(a),
        'force_rmse_vector_norm_eV_A':float(np.sqrt(np.mean(np.sum(f*f,axis=1)))),
        'max_force_norm_eV_A':float(np.linalg.norm(f,axis=1).max()),
        'stress_eV_A3':stress.tolist(),'scf_converged':conv,'finite_energy_force_stress':finite,
        'scf_iterations':int(a.calc.get_number_of_iterations()),
        'run_status':'complete_candidate_only' if conv and finite else 'rejected_unconverged_or_nonfinite'})
    if conv and finite:
        a.info['energy_force_consistent_free_eV']=e_free
        a.info['energy_extrapolated_eV']=e_extrap
        a.info['v14_candidate_id']=m['candidate_id']
        a.info['dft_label_status']='candidate_only_not_yet_accepted_for_fit'
        a.calc=SinglePointCalculator(a,energy=e_free,forces=f,stress=stress.reshape(-1))
        write(label,a,format='extxyz')
except Exception as e:
    rec.update({'run_status':'failed','error':f'{type(e).__name__}: {e}','traceback':traceback.format_exc()})
rec['elapsed_seconds']=round(time.time()-t0,3)
for name,path in [('geometry_sha256',geom),('log_sha256',log),('label_sha256',label)]:
    if path.exists():rec[name]=hashlib.sha256(path.read_bytes()).hexdigest()
result.write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps({k:rec.get(k) for k in ('candidate_id','run_status','scf_iterations','elapsed_seconds','energy_fit_target_eV','max_force_norm_eV_A')},indent=2),flush=True)
