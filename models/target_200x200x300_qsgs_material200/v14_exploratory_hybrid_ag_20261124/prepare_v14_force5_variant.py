#!/workspace/venvs/gpaw/bin/python
"""Create an isolated V14 force-emphasis sensitivity fit."""
import hashlib, json, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BASE = ROOT / 'fit_v14_targeted_ag'
OUT = ROOT / 'fit_v14_targeted_ag_force5'
if OUT.exists():
    raise SystemExit(f'Refusing to overwrite {OUT}')
shutil.copytree(BASE / 'training', OUT / 'training')
shutil.copytree(BASE / 'validation', OUT / 'validation')
shutil.copy2(BASE / 'AgTiSiC_POD_v14.pod', OUT / 'AgTiSiC_POD_v14.pod')
(OUT / 'fitpod.in').write_text('clear\nunits metal\natom_style atomic\nfitpod AgTiSiC_POD_v14.pod fitpod.data\n')
(OUT / 'fitpod.data').write_text(f'''file_format extxyz
file_extension xyz
path_to_training_data_set {OUT / 'training'}
path_to_test_data_set {OUT / 'validation'}
basename_for_output_files v14_targeted_ag_f5
fitting_weight_energy 100.0
fitting_weight_force 5.0
fitting_weight_stress 0.0
error_analysis_for_training_data_set 1
error_analysis_for_test_data_set 1
energy_force_calculation_for_training_data_set 0
energy_force_calculation_for_test_data_set 0
fraction_training_data_set 1.0
fraction_test_data_set 1.0
randomize_training_data_set 0
randomize_test_data_set 0
fitting_regularization_parameter 0.1
save_pod_descriptors 0
compute_pod_descriptors 0
''')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
record = {'status': 'sensitivity variant staged', 'parent_fit_manifest_sha256': sha(BASE/'fit_data_manifest.json'),
          'training_count': 31, 'validation_count': 6,
          'fitting_weight_energy': 100.0, 'fitting_weight_force': 5.0,
          'note': 'Only force fitting weight differs from V14 targeted baseline; no physical data or descriptor changes.'}
(OUT / 'variant_manifest.json').write_text(json.dumps(record, indent=2)+'\n')
print(json.dumps({'fit_dir': str(OUT), **record}, indent=2))
