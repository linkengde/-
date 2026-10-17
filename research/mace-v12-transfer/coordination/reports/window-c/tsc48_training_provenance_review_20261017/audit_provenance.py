"""Read-only audit of existing Git blobs; never executes the DFT runner."""
import io,json,hashlib,subprocess
from pathlib import Path
import numpy as np
from ase.io import read
C='fa699e1340e8777a7095bb8c639a7247c28d1956'
R='research/mace-v12-transfer/coordination/reports/window-e/dft_supplement_results_20261016/'
S=R+'tscbop_run0_sixpoint/analysis/training_use_provenance_supplement_20261017/'
O=Path(__file__).resolve().parent
sources={}
def blob(p):
 b=subprocess.check_output(['git','show',C+':'+p]);sources[p]=hashlib.sha256(b).hexdigest();return b
def js(p):return json.loads(blob(p))
lines=blob(S+'SHA256SUMS.txt').decode().splitlines();count=0
for line in lines:
 if not line.strip() or line.startswith('#'):continue
 h,p=line.split(maxsplit=1);assert hashlib.sha256(blob(p)).hexdigest()==h;count+=1
assert count==16
sup=js(S+'E_TO_C_TRAINING_PROVENANCE_SUPPLEMENT.json');lock=js(R+'EXECUTION_LOCK_K6X6X2_MSR1.json');runner=blob(R+'run_single_point_msr1.py');pre=js(R+'PRELAUNCH_PREFLIGHT_MSR1_K6X6X2_COLLECTORFIX.json');paw=js(S+'PAW_SETUP_SHA256.json')
lockh=sources[R+'EXECUTION_LOCK_K6X6X2_MSR1.json'];runh=hashlib.sha256(runner).hexdigest()
assert lockh==sup['execution_lock']['sha256']==pre['parameter_lock_sha256']
assert runh==sup['runner']['sha256']==pre['runner_sha256']
cm=js('research/mace-v12-transfer/coordination/reports/window-e/dft_supplement_prep_20261015/ti3sic2_bulk_p0/structures/DFT_INPUT_CANDIDATE_MANIFEST.json')
rt=js('research/mace-v12-transfer/coordination/reports/window-e/dft_supplement_prep_20261015/E_RUNTIME_PREFLIGHT.json')
for elem in ['Ti','Si','C']:
 n=elem+'.PBE.gz';h=pre['PAW_sha256'][n]
 assert all(v==h for k,v in paw['files'][n].items() if k.endswith('_sha256'))
 assert rt['PAW_setup_sha256'][n]==cm['proposed_calculation']['PAW_SHA256'][n]==h
for rec in sup['six_commands']:
 cmd=js(rec['path']);args=cmd['argv'];get=lambda k:args[args.index(k)+1]
 assert cmd['runner_sha256']==runh and cmd['parameter_lock_sha256']==lockh
 assert get('--execution-lock-sha256')==lockh and get('--input-sha256')==rec['input_sha256']
 assert get('-n')=='4' and cmd['omp_num_threads']=='1'
 inp=str(Path(rec['path']).parent/'input.extxyz');assert hashlib.sha256(blob(inp)).hexdigest()==rec['input_sha256']
assert lock['scientific_method']['kpoints']['size']==[6,6,2]
assert b"kpts=tuple(method['kpoints']['size'])" in runner
assert b"atoms.get_potential_energy(force_consistent=False)" in runner and b"atoms.get_potential_energy(force_consistent=True)" in runner
m=js(sup['mesh_comparison']['path']);base='research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/calculations/'
arr=[];summaries=[]
for grid in ['k6x6x2','k6x6x4']:
 name='Ti3SiC2_baseline_'+grid+'_v19';p=base+name+'/'+name+'_PW_PBE.extxyz';b=blob(p)
 assert hashlib.sha256(b).hexdigest()==m['archive_extxyz_sha256'][grid]
 arr.append(read(io.StringIO(b.decode()),format='extxyz'));summaries.append(js(base+name+'/summary.json'))
a,b=arr;assert np.array_equal(a.numbers,b.numbers) and np.allclose(a.positions,b.positions,atol=1e-8,rtol=0) and np.allclose(a.cell,b.cell,atol=1e-8,rtol=0)
df=b.arrays['PW_PBE_forces']-a.arrays['PW_PBE_forces']
measured={'native_delta_meV_atom':(b.info['PW_PBE_energy_eV']-a.info['PW_PBE_energy_eV'])*1000/48,'free_energy_delta_meV_atom':(summaries[1]['free_energy_eV_cell']-summaries[0]['free_energy_eV_cell'])*1000/48,'force_vector_rms_eV_A':float(np.sqrt(np.mean(np.sum(df*df,axis=1)))),'force_vector_max_eV_A':float(np.linalg.norm(df,axis=1).max())}
for k,v in measured.items():assert abs(v-sup['mesh_comparison'][k])<1e-8
result={'date':'2026-10-17','E_COMMIT':C,'SUPPLEMENT_SHA_CHECKS':count,'LOCK_RUNNER_BINDING':'PASS','PAW_ARCHIVE_CONSISTENCY':'PASS_WITH_BATCH_LEVEL_LIMITATION','ACTUAL_KPOINTS':[6,6,2],'KPOINT_METADATA_CONFLICT':'RESOLVED_AS_STALE_NONEXECUTABLE_TEXT','PARENT_MESH_RECOMPUTATION':measured,'PERTURBED_GEOMETRY_CONVERGENCE':'NOT_DEMONSTRATED','STRESS_CONVERGENCE':'NOT_DEMONSTRATED','PROVENANCE_REVIEW':'ACCEPTED_WITH_LIMITATIONS','TRAIN_USE_APPLICATION_READY':True,'TRAIN_USE':'HOLD_PENDING_EXPLICIT_USE_APPROVAL_AND_APPROPRIATE_NUMERICAL_SCOPE','TRAIN':[],'DEV':[],'TEST':[],'role_change_authorized':False,'fit_authorized':False,'new_calculation_authorized':False,'EXECUTION_LOCK_SHA256':lockh,'RUNNER_SHA256':runh,'PAW_evidence_scope':'archived batch preflight plus matching archived source and later E attestation; not C live or per-job/per-rank hashing'}
(O/'REVIEW_DECISION.json').write_text(json.dumps(result,indent=2)+'\n')
(O/'SOURCE_SHA256.json').write_text(json.dumps({'commit':C,'files':sources},indent=2)+'\n')
print(json.dumps(result,indent=2))
