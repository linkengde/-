import csv
import json
import pathlib
import subprocess
import sys
import tempfile

repo = pathlib.Path('/workspace/-')
report_dir = repo / 'research/mace-v12-transfer/coordination/reports/window-a/v19_Ti3SiC2_bulk_mesh_comparison'
pair_script = report_dir / 'compare_mesh_pair.py'
pairs = [
    ('Ti3SiC2_baseline_k4x4x2_v19', 'Ti3SiC2_baseline_k6x6x2_v19'),
    ('Ti3SiC2_baseline_k6x6x2_v19', 'Ti3SiC2_baseline_k8x8x2_v19'),
    ('Ti3SiC2_baseline_k8x8x2_v19', 'Ti3SiC2_baseline_k10x10x2_v19'),
    ('Ti3SiC2_baseline_k6x6x2_v19', 'Ti3SiC2_baseline_k10x10x2_v19'),
    ('Ti3SiC2_baseline_k6x6x4_v19', 'Ti3SiC2_baseline_k8x8x4_v19'),
    ('Ti3SiC2_baseline_k8x8x4_v19', 'Ti3SiC2_baseline_k10x10x4_v19'),
    ('Ti3SiC2_baseline_k6x6x4_v19', 'Ti3SiC2_baseline_k10x10x4_v19'),
]
rows = []
for left, right in pairs:
    with tempfile.NamedTemporaryFile(suffix='.json') as tmp:
        subprocess.run([
            sys.executable, str(pair_script), '--left', left, '--right', right,
            '--output', tmp.name,
        ], cwd=repo, check=True, stdout=subprocess.DEVNULL)
        result = json.loads(pathlib.Path(tmp.name).read_text())
    row = {
        'left': left,
        'right': right,
        'source_hash': result['source_hash'],
        'same_geometry_method_except_kpts': result['same_geometry_method_except_kpts'],
        'archive_verification': result['archive_verification'],
        'native_energy_delta_meV_atom': result['native_energy_delta_meV_atom'],
        'free_energy_delta_meV_atom': result['free_energy_delta_meV_atom'],
        'force_vector_rms_eV_A': result['force']['vector_rms_eV_A'],
        'maximum_atom_vector_difference_eV_A': result['force']['maximum_atom_vector_difference_eV_A'],
        'pass': result['pass'],
    }
    if not row['same_geometry_method_except_kpts']:
        raise SystemExit(f'Unexpected non-matched pair: {left} vs {right}')
    if result['archive_verification']['left'] != 'PASS' or result['archive_verification']['right'] != 'PASS':
        raise SystemExit(f'Unverified archive in pair: {left} vs {right}')
    rows.append(row)
report = {
    'status': 'EXISTING_FIXED_GEOMETRY_MESH_MATRIX',
    'scope': 'Same 48-atom frozen Ti3SiC2 baseline and PBE/PW500 settings; only k-point meshes vary. Not interface validation or thermal stability.',
    'thresholds': {
        'absolute_native_energy_delta_meV_atom': 2.0,
        'absolute_free_energy_delta_meV_atom': 2.0,
        'force_vector_rms_eV_A': 0.01,
    },
    'pairs': rows,
    'interpretation': [
        'All archived pair endpoints passed SHA and archive verification.',
        'At kz=2, 8x8 to 10x10 and cumulative 6x6 to 10x10 pass all budgets; 4x4 to 6x6 and 6x6 to 8x8 exceed the force RMS budget.',
        'At kz=4, 8x8 to 10x10 and cumulative 6x6 to 10x10 pass all budgets; 6x6 to 8x8 narrowly exceeds the force RMS budget.',
        'No further pristine-baseline mesh point is justified by these comparisons alone; perturbed/interface geometries require their own force checks.',
    ],
}
json_path = report_dir / 'mesh_matrix_kz2_kz4.json'
json_path.write_text(json.dumps(report, indent=2) + '\n')
csv_path = report_dir / 'mesh_matrix_kz2_kz4.csv'
with csv_path.open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=[
        'left', 'right', 'native_energy_delta_meV_atom', 'free_energy_delta_meV_atom',
        'force_vector_rms_eV_A', 'maximum_atom_vector_difference_eV_A',
        'native_energy_pass', 'free_energy_pass', 'force_vector_rms_pass',
    ])
    writer.writeheader()
    for row in rows:
        writer.writerow({
            'left': row['left'], 'right': row['right'],
            'native_energy_delta_meV_atom': row['native_energy_delta_meV_atom'],
            'free_energy_delta_meV_atom': row['free_energy_delta_meV_atom'],
            'force_vector_rms_eV_A': row['force_vector_rms_eV_A'],
            'maximum_atom_vector_difference_eV_A': row['maximum_atom_vector_difference_eV_A'],
            'native_energy_pass': row['pass']['native_energy'],
            'free_energy_pass': row['pass']['free_energy'],
            'force_vector_rms_pass': row['pass']['force_vector_rms'],
        })
print(json.dumps({'pairs': len(rows), 'outputs': [str(json_path), str(csv_path)], 'thresholds': report['thresholds']}, indent=2))
