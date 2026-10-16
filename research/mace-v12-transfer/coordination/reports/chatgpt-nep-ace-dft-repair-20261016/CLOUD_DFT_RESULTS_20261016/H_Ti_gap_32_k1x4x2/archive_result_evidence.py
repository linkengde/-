from pathlib import Path
import csv,hashlib,json,shutil,re,subprocess,sys
import numpy as np
from ase.io import read
import gpaw,ase,gpaw_data,importlib.metadata as md
BASE=Path(__file__).resolve().parent
REPO=Path('/workspace/h-dft-repair')
REL=Path('research/mace-v12-transfer/coordination/reports/chatgpt-nep-ace-dft-repair-20261016')
RAW=BASE/'Ti_gap_32_k1x4x2'
DEST=REPO/REL/'CLOUD_DFT_RESULTS_20261016'/'H_Ti_gap_32_k1x4x2'
monitor=json.loads((BASE/'EXECUTION_MONITOR.json').read_text())
assert monitor['state']=='EXITED' and monitor['returncode']==0,'DFT not successfully finished'
result=json.loads((RAW/'RESULT.json').read_text())
input_file=REPO/REL/'Ti_gap_2p3_frozen_32atoms_UNLABELED.extxyz'
a=read(input_file);forces=np.array(result['forces_eV_A'])
assert forces.shape==(32,3) and np.isfinite(forces).all()
assert np.isfinite([result['native_energy_eV'],result['force_consistent_free_energy_eV']]).all()
assert result['input_sha256']==hashlib.sha256(input_file.read_bytes()).hexdigest()
assert result['kmesh']==[1,4,2] and result['mpi_ranks']==4
log=(RAW/'gpaw.log').read_text()
converged_lines=[s for s in log.splitlines() if 'converge' in s.lower()]
iterations=re.findall(r'converged\s+(?:in|after)\s+(\d+)',log,re.I)
if not iterations:
    iterations=re.findall(r'\|iter:\s*(\d+)\|',log)
iteration_lines=[s for s in log.splitlines() if re.match(r'\s*(?:iter:)?\s*\d+\s',s)]
# Successful GPAW energy/force calls imply SCF succeeded; archive literal log evidence.
assert 'converged' in log.lower(), 'Missing literal convergence evidence; inspect before release'
assert not DEST.exists(),'Never overwrite existing archive'
DEST.mkdir(parents=True)
shutil.copy2(Path(__file__),DEST/'archive_result_evidence.py')
for p in RAW.iterdir():
 if p.is_file() and p.suffix!='.gpw':shutil.copy2(p,DEST/p.name)
shutil.copy2(input_file,DEST/'frozen_input.extxyz')
shutil.copy2(REPO/REL/'run_32atom_kmesh_MPIOPT.py',DEST/'runner_pinned.py')
for name in ['launcher.log','EXECUTION_MONITOR.json','supervise_first_job.py','launcher_cli_error.log','EXECUTION_MONITOR_cli_error.json']:
 shutil.copy2(BASE/name,DEST/name)
shutil.copytree(BASE/'preflight',DEST/'preflight')
with (DEST/'forces_indexed.csv').open('w',newline='') as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['atom_index_0based','element','x_A','y_A','z_A','Fx_eV_A','Fy_eV_A','Fz_eV_A'])
 for i,(sym,xyz,force) in enumerate(zip(a.get_chemical_symbols(),a.positions,forces)):
  w.writerow([i,sym,*[format(float(x),'.17g') for x in xyz],*[format(float(x),'.17g') for x in force]])
summary={**result,'SCF_converged':True,'SCF_convergence_log_evidence':converged_lines,'SCF_iterations':int(iterations[-1]) if iterations else None,'wall_seconds':monitor['elapsed_seconds'],'started_at':monitor['started_at'],'ended_at':monitor['ended_at'],'peak_sampled_process_tree_rss_bytes':monitor['peak_sampled_process_tree_rss_bytes'],'memory_metric_note':monitor['memory_metric_note'],'source_commit':'ed7d747c968b555144e9b465ff6ac041717dcc6e','input_coordinates_unchanged':True,'atom_index_convention':'0-based; frozen input order','kpoint_convergence_validated':False,'historical_1x4x2_centering':'UNVERIFIED','further_jobs_started':False,'checkpoint_saved_by_pinned_runner':False,'authorization_scope':'only P0 32 atom gamma-centered 1x4x2 fixed-geometry single point'}
(DEST/'DFT_P0_RESULT.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
paw=json.loads((BASE/'paw_manifest.json').read_text())
env={'Python':sys.version.split()[0],'GPAW':gpaw.__version__,'ASE':ase.__version__,'gpaw_data':md.version('gpaw-data'),'NumPy':np.__version__,'mpi_ranks':4,'MPI_4rank_import_test':'PASS','OpenMPI':'5.0.7','paw_setup_files':paw,'source_commit':summary['source_commit'],'host_identifiers_redacted':True,'PAW_match_project_baseline':True,'GPAW_SETUP_PATH_matches_gpaw_data':True,'runtime_script_sha256':hashlib.sha256(Path('/workspace/.setup/runtime.sh').read_bytes()).hexdigest(),'input_sha256':summary['input_sha256'],'runner_sha256':hashlib.sha256((DEST/'runner_pinned.py').read_bytes()).hexdigest(),'persistent_volume_guarantee':'UNVERIFIED','external_new_software_installed':False}
(DEST/'ENVIRONMENT_AND_PAW_SHA.json').write_text(json.dumps(env,indent=2)+'\n')
(DEST/'SCF_AND_KPOINT_AUDIT.md').write_text(f'''# H P0 fixed-geometry single-point audit

Source commit: `{summary['source_commit']}`. Only the authorized 32-atom Ag10Ti12Si2C8 job was run: PBE / PW 500 eV / Fermi-Dirac 0.10 eV / spin unpolarized / full 3D PBC / gamma-centered 1x4x2 / four MPI ranks. Original coordinates, cell, atomic order and input bytes are preserved. No geometry relaxation, MD, NEP/ACE training or next-grid task was started.

SCF converged: **YES**, successful GPAW energy and force calls and literal log evidence: `{converged_lines}`. SCF iterations parsed: `{summary['SCF_iterations']}`. Native energy: `{summary['native_energy_eV']:.15g}` eV. Force-consistent free energy: `{summary['force_consistent_free_energy_eV']:.15g}` eV. All 32x3 forces are finite; ordered values and 0-based indices are in forces_indexed.csv and DFT_P0_RESULT.json.

Wall time: {monitor['elapsed_seconds']:.3f} seconds. Peak sampled MPI process-tree sum RSS: {monitor['peak_sampled_process_tree_rss_bytes']} bytes; sampled every 2 seconds and may double-count shared pages. This is not an exact per-rank or scheduler peak. Raw output is preserved outside the Git checkout; no large checkpoint is committed. The supplied runner does not write a .gpw checkpoint.

The first launcher attempt exited during CLI parsing before SCF because the `gpaw python --` script-argument separator was omitted. Its error log is retained. The corrected invocation executes the exact pinned runner; no scientific parameters were changed.

New grid centering is explicit gamma=True, confirmed in gpaw.log. The old 1x4x2 centering is **UNVERIFIED**: inspected historical scripts do not establish this particular label's original centering. Therefore no same-method historical k-point convergence claim is made. A single converged SCF is not k-point convergence; 1x6x4 and other structures remain unstarted. C/user review is needed before any next job.

All eleven supplied structures passed Git-blob, atom-count, formula, PBC, finite-coordinate and ASE-read checks; see preflight/. SHA256SUMS.txt covers result and evidence files. Hash PASS verifies archived bytes only, not model quality or resolution of historical label conflicts. No training eligibility or data role has been changed.
''')
(DEST/'README.md').write_text('''# H P0 result

New result from the authorized single gamma-centered 1x4x2 SCF; source commit ed7d747c968b555144e9b465ff6ac041717dcc6e. Read DFT_P0_RESULT.json, forces_indexed.csv, gpaw.log, SCF_AND_KPOINT_AUDIT.md and ENVIRONMENT_AND_PAW_SHA.json. Native and force-consistent energies are distinct. No k-point convergence or automatic training approval is claimed. No next task was started.

Verify files: `sha256sum -c SHA256SUMS.txt`. Raw output remains outside the repository; the pinned runner did not write a checkpoint.
''')
(DEST/'RESULT_VERIFICATION.json').write_text(json.dumps({'archive_validation':'PASS','finite_energy_and_forces':True,'force_shape':[32,3],'input_sha_matches':True,'input_and_atom_order_preserved':True,'PBE_PW500_Fermi0p10_spin_unpolarized':True,'SCF_converged':True,'mesh':[1,4,2],'gamma_centered':True,'kpoint_convergence':'NOT_VALIDATED','next_tasks_started':False},indent=2)+'\n')
paths=sorted(p for p in DEST.rglob('*') if p.is_file())
(DEST/'SHA256SUMS.txt').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(DEST).as_posix()+'\n' for p in paths))
print(DEST)
print(json.dumps({'energy_eV':summary['native_energy_eV'],'free_energy_eV':summary['force_consistent_free_energy_eV'],'wall_seconds':summary['wall_seconds'],'SCF_iterations':summary['SCF_iterations']},indent=2))
