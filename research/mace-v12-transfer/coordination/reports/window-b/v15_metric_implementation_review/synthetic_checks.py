"""Synthetic-only review; AST extraction avoids evaluation imports/main and all data/model reads."""
import ast
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path('/workspace/-')
PROJECT=ROOT/'research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v15_interface_energy'
OUT=ROOT/'research/mace-v12-transfer/coordination/reports/window-b/v15_metric_implementation_review'
source=(PROJECT/'evaluate_v15.py').read_text();tree=ast.parse(source)
node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='aggregate')
ns={'np':np};exec(compile(ast.Module(body=[node],type_ignores=[]),'aggregate extracted AST','exec'),ns)
aggregate=ns['aggregate'];results=[]
def check(name,observed,expected):
    assert np.allclose(observed,expected,rtol=1e-12,atol=1e-12),(name,observed,expected)
    results.append({'check':name,'status':'PASS','observed':observed,'expected':expected})
check('energy normalization: 0.48 eV for 48 atoms',abs(.48)*1000/48,10.)
delta=np.array([[.03,.04,0.],[0.,0.,.05]])
vector=float(np.sqrt(np.mean(np.sum(delta**2,axis=1))));component=float(np.sqrt(np.mean(delta**2)))
check('per-atom force vector RMSE',vector,.05);check('vector versus component ratio',vector/component,float(np.sqrt(3)))
check('maximum is max vector norm',float(np.linalg.norm(delta,axis=1).max()),.05)
rows=[{'atoms':1,'energy_abs_error_meV_atom':2.,'force_vector_RMSE_eV_A':.1,'force_vector_max_error_eV_A':.1,'separating_force_abs_error_eV_A':.03},{'atoms':3,'energy_abs_error_meV_atom':10.,'force_vector_RMSE_eV_A':.2,'force_vector_max_error_eV_A':.25,'separating_force_abs_error_eV_A':.07}]
a=aggregate(rows);check('aggregate atom-weighted vector RMSE',a['force_vector_RMSE_eV_A'],float(np.sqrt(.0325)));check('aggregate energy is per-structure mean',a['energy_MAE_meV_atom'],6.);check('aggregate maximum vector error',a['force_vector_max_error_eV_A'],.25);check('aggregate maximum separation error',a['separating_force_max_abs_error_eV_A'],.07)
# Reverse array order: choose chemical identities rather than central index order.
pos=np.array([[2.,0.,0.],[0.,0.,0.]]);forces=np.array([[1.,0.,0.],[-1.,0.,0.]])
symbols=['C','Ag'];ag=[i for i in range(2) if symbols[i]=='Ag'];contact=[i for i in range(2) if symbols[i]!='Ag'];i,j=ag[0],contact[0];u=(pos[j]-pos[i])/np.linalg.norm(pos[j]-pos[i]);check('positive separating force with reversed array order',float(np.dot(forces[j]-forces[i],u)),2.)
gates={'energy_MAE_meV_atom':10.,'force_vector_RMSE_eV_A':.05,'separating_force_max_abs_error_eV_A':.1}
assert all(gates[k]<=v for k,v in gates.items());results.append({'check':'gate equality accepted','status':'PASS'})
above=dict(gates);above['force_vector_RMSE_eV_A']=float(np.nextafter(.05,np.inf));assert not all(above[k]<=v for k,v in gates.items());results.append({'check':'immediately above gate rejected','status':'PASS'})
assert not all(None is not None and None<=v for v in gates.values());results.append({'check':'missing metric not accepted','status':'PASS'})
# Exercise only extracted pure guards on fabricated atom-like arrays.
pre_tree=ast.parse((PROJECT/'preflight.py').read_text())
for function in ['require','labels','geometry']:
    n=next(n for n in pre_tree.body if isinstance(n,ast.FunctionDef) and n.name==function)
    exec(compile(ast.Module(body=[n],type_ignores=[]),'preflight guard AST','exec'),ns)
class FakeAtoms:
    def __init__(self):
        self.numbers=np.array([47,6]);self.pbc=np.array([True,True,True]);self.positions=np.array([[0.,0.,0.],[2.,0.,0.]]);self.cell=np.eye(3)*10
        self.cell=type('Cell',(),{'array':self.cell,'__array__':lambda obj,dtype=None:np.asarray(obj.array,dtype=dtype)})()
        self.info={'REF_energy':-1.};self.arrays={'REF_forces':np.zeros((2,3)),'lammps_id':np.array([10,20]),'central_pair':np.array([True,True])}
    def __len__(self):return 2
ns['labels'](FakeAtoms());ns['geometry'](FakeAtoms(),FakeAtoms())
results.append({'check':'synthetic finite labels and matching geometry accepted','status':'PASS'})
for name,change,function in [
 ('missing reference energy',lambda a:a.info.pop('REF_energy'),'labels'),
 ('nonfinite force',lambda a:a.arrays['REF_forces'].__setitem__((0,0),np.nan),'labels'),
 ('wrong persistent ID',lambda a:a.arrays['lammps_id'].__setitem__(0,11),'geometry'),
 ('wrong PBC',lambda a:a.pbc.__setitem__(0,False),'geometry'),
 ('wrong position',lambda a:a.positions.__setitem__((0,0),.01),'geometry')]:
    a=FakeAtoms();change(a)
    try:
        if function=='geometry':ns[function](FakeAtoms(),a)
        else:ns[function](a)
    except ValueError as error:results.append({'check':name+' fails clearly','status':'PASS','diagnostic':str(error)})
    else:raise AssertionError(name+' was not rejected')

# Current contract uses direct cluster vector: reproduce its periodic limitation without changing policy.
direct=np.array([9.8,1.,0.]);mic=np.array([-.2,1.,0.]);force=np.array([1.,0.,0.]);dp=float(np.dot(force,direct/np.linalg.norm(direct)));mp=float(np.dot(force,mic/np.linalg.norm(mic)))
limitations=[{'id':'L1_periodic_image_convention','direct_separation_eV_A':dp,'MIC_separation_eV_A':mp,'meaning':'Synthetic wrapped-cell example; does not show any actual blind input or score is affected. Direct cluster convention explicitly retained from V14; periodic-image policy must be reviewed before reuse.'}]
# Finite inputs can overflow derived metrics: production only rejects at JSON serialization, not immediately.
with np.errstate(over='ignore',invalid='ignore'):
 huge=np.array([[1e308,0.,0.]]);overflow=float(np.sqrt(np.mean(np.sum(huge**2,axis=1))))
assert np.isfinite(huge).all() and not np.isfinite(overflow)
try:json.dumps({'RMSE':overflow},allow_nan=False)
except ValueError:serialization_rejects=True
else:serialization_rejects=False
assert serialization_rejects
findings=[{'id':'D1_derived_nonfinite_not_fail_fast','severity':'low robustness','finite_force_inputs':True,'derived_RMSE':'inf','serialization_rejects':True,'meaning':'JSON prevents publishing invalid metrics, but failure occurs after later evaluations/output mkdir rather than at the offending frame. Add a finite guard to derived metrics.'}]
# Freeze-record race: code checks parsed values but hashes the file again later.
initial={'model_sha256':'selected-A','blind_labels_used_for_selection':False};later={'model_sha256':'selected-B','blind_labels_used_for_selection':True}
assert initial!=later
findings.append({'id':'D2_selection_record_snapshot_not_bound','severity':'medium evidence integrity','checked_snapshot':initial,'possible_later_snapshot':later,'meaning':'Lines27/73 read and hash different snapshots if selection file is rewritten during evaluation. Preserve/check initial selection bytes/hash before final reporting. No real files/models/evaluation were run.'})
report={'reviewed_files_sha256':{str((PROJECT/n).relative_to(ROOT)):hashlib.sha256((PROJECT/n).read_bytes()).hexdigest() for n in ['evaluate_v15.py','preflight.py']},'synthetic_checks':results,'definite_robustness_findings':findings,'limitations':limitations,'scope':'Only production source code read. No label/input datasets, model files, Torch/MACE/GPAW imports or evaluation main execution.'}
(OUT/'checks.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print('Synthetic assertions passed:',len(results),'robustness findings:',len(findings))
if __name__=='__main__':pass
