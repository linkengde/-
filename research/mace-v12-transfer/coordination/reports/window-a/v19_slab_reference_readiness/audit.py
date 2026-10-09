#!/usr/bin/env python3
"""Read-only verification of existing AgSi slab numerical controls."""
import hashlib
import json
from pathlib import Path
import numpy as np
from ase.io import read
from ase.geometry import find_mic

ROOT = Path('/workspace/-')
BASE = ROOT / 'research/mace-v12-transfer/periodic_interface_v4'
REPORT = ROOT / 'research/mace-v12-transfer/coordination/reports/window-a/v19_slab_reference_readiness'

ARCHIVES = {
    'AgSi_k6_sigma010': (BASE / 'pbe_interface_v19_cif_convergence_pilot', 'AgSi_COD9009647_pilot_k6x6'),
    'AgSi_k8_sigma010': (BASE / 'pbe_interface_v19_cif_evenmesh_pilot', 'AgSi_COD9009647_pilot_k8x8'),
    'AgSi_vacuum26_k6_sigma010': (BASE / 'pbe_interface_v19_cif_vacuum_pilot', 'AgSi_COD9009647_vacuum26_k6x6_sigma0p10'),
    'dipole_off_16A': (BASE / 'pbe_interface_v19_cif_slab_controls', 'AgSi_originalvac_slab_dipole_off_v19'),
    'dipole_on_16A': (BASE / 'pbe_interface_v19_cif_slab_controls', 'AgSi_originalvac_slab_dipole_on_v19'),
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def verify_archive(key):
    folder, label = ARCHIVES[key]
    folder = folder / 'calculations' / label
    checks = {}
    for line in (folder / 'SHA256SUMS.txt').read_text().splitlines():
        digest, name = line.split('  ', 1)
        checks[name] = sha(folder / name) == digest
    verification = json.loads((folder / 'verification.json').read_text())
    summary = json.loads((folder / 'summary.json').read_text())
    record = next(r for r in json.loads((ARCHIVES[key][0] / 'input_manifest.json').read_text())['records'] if r['label'] == label)
    atoms = read(folder / f'{label}_PW_PBE.extxyz', format='extxyz')
    if not all(checks.values()) or verification.get('status') != 'PASS' or not all(verification.get('checks', {}).values()):
        raise SystemExit(f'Archive integrity/verifier failed: {label}')
    if summary.get('scf_converged') is not True or summary.get('source_sha256') != record['input_sha256']:
        raise SystemExit(f'Archive convergence/source mismatch: {label}')
    vacuum = float(atoms.cell.lengths()[2] - np.ptp(atoms.positions[:, 2]))
    return {
        'label': label,
        'archive_status': verification['status'],
        'archive_sha_checks': checks,
        'input_sha256': record['input_sha256'],
        'kpts': list(summary['method']['kpts']),
        'smearing_eV': summary['method']['smearing_eV'],
        'pbc': atoms.pbc.tolist(),
        'poissonsolver': summary['method'].get('poissonsolver', {}),
        'cell_z_A': float(atoms.cell.lengths()[2]),
        'coordinate_span_z_A': float(np.ptp(atoms.positions[:, 2])),
        'geometric_empty_length_A': vacuum,
        'native_energy_eV': summary['energy_eV_cell'],
        'free_energy_eV': summary['free_energy_eV_cell'],
        'iterations': summary['scf_iterations'],
        'max_force_eV_A': summary['fmax_eV_A'],
        'summary': summary,
        'atoms': atoms,
    }

def pair_metrics(left, right, pair_projection=True):
    a, sa = left['atoms'], left['summary']
    b, sb = right['atoms'], right['summary']
    equal = {
        'symbols': np.array_equal(a.numbers, b.numbers),
        'positions': np.array_equal(a.positions, b.positions),
        'cell': np.array_equal(a.cell.array, b.cell.array),
        'atom_ids': np.array_equal(a.arrays['lammps_id'], b.arrays['lammps_id']),
        'pair_marker': np.array_equal(a.arrays['central_pair'], b.arrays['central_pair']),
    }
    if not all(equal.values()):
        raise SystemExit(f'Pair geometry/IDs mismatch: {equal}')
    f0, f1 = a.arrays['PW_PBE_forces'], b.arrays['PW_PBE_forces']
    df = f1 - f0
    vectors = np.linalg.norm(df, axis=1)
    n = len(a)
    result = {
        'same_geometry_ids': equal,
        'native_energy_delta_meV_atom': float((sb['energy_eV_cell'] - sa['energy_eV_cell']) * 1000 / n),
        'free_energy_delta_meV_atom': float((sb['free_energy_eV_cell'] - sa['free_energy_eV_cell']) * 1000 / n),
        'force_vector_rms_eV_A': float(np.sqrt(np.mean(vectors ** 2))),
        'maximum_atom_vector_difference_eV_A': float(vectors.max()),
        'marked_pair_projection_delta_eV_A': None,
        'budgets': {'absolute_energy_meV_atom': 2.0, 'force_vector_rms_eV_A': 0.01, 'absolute_pair_projection_eV_A': 0.02},
    }
    if pair_projection:
        marked = a.arrays['central_pair'].astype(bool)
        symbols = np.asarray(a.get_chemical_symbols())
        i, j = np.flatnonzero(marked)
        if a[i].symbol != 'Ag':
            i, j = j, i
        temp = a.copy(); temp.set_pbc([True, True, False])
        direction = temp.get_distance(int(i), int(j), mic=True, vector=True)
        direction /= np.linalg.norm(direction)
        result['marked_pair_projection_delta_eV_A'] = float(np.dot(df[j] - df[i], direction))
    result['pass'] = {
        'native_energy': abs(result['native_energy_delta_meV_atom']) <= 2.0,
        'free_energy': abs(result['free_energy_delta_meV_atom']) <= 2.0,
        'force_vector_rms': result['force_vector_rms_eV_A'] <= 0.01,
        'marked_pair_projection': result['marked_pair_projection_delta_eV_A'] is None or abs(result['marked_pair_projection_delta_eV_A']) <= 0.02,
    }
    return result

def vacuum_pair_metrics(left, right):
    """Compare the 16 -> 26 A control with a rigid +5 A z translation."""
    a, sa = left['atoms'], left['summary']
    b, sb = right['atoms'], right['summary']
    shift = b.positions - a.positions
    equal = {
        'symbols': np.array_equal(a.numbers, b.numbers),
        'rigid_global_translation_plus_5A_z': bool(np.allclose(shift, [0.0, 0.0, 5.0], atol=1e-10, rtol=0.0)),
        'xy_cell': bool(np.allclose(a.cell.array[:2], b.cell.array[:2], atol=1e-12, rtol=0.0)),
        'z_cell_plus_10A': abs(float(b.cell.lengths()[2] - a.cell.lengths()[2]) - 10.0) <= 1e-10,
        'pbc': np.array_equal(a.pbc, b.pbc),
        'atom_ids': np.array_equal(a.arrays['lammps_id'], b.arrays['lammps_id']),
        'pair_marker': np.array_equal(a.arrays['central_pair'], b.arrays['central_pair']),
    }
    if not all(equal.values()):
        raise SystemExit(f'Vacuum pair geometry/IDs mismatch: {equal}')
    if sa['method'] != sb['method']:
        raise SystemExit('Vacuum pair method mismatch')
    f0, f1 = a.arrays['PW_PBE_forces'], b.arrays['PW_PBE_forces']
    df = f1 - f0
    vectors = np.linalg.norm(df, axis=1)
    i, j = np.flatnonzero(a.arrays['central_pair'].astype(bool))
    if a[i].symbol != 'Ag':
        i, j = j, i
    temp = a.copy(); temp.set_pbc([True, True, False])
    direction = temp.get_distance(int(i), int(j), mic=True, vector=True)
    direction /= np.linalg.norm(direction)
    result = {
        'same_method': True,
        'same_method_except_vacuum_geometry': True,
        'geometry_transform_checks': equal,
        'native_energy_delta_meV_atom': float((sb['energy_eV_cell'] - sa['energy_eV_cell']) * 1000 / len(a)),
        'free_energy_delta_meV_atom': float((sb['free_energy_eV_cell'] - sa['free_energy_eV_cell']) * 1000 / len(a)),
        'force_vector_rms_eV_A': float(np.sqrt(np.mean(vectors ** 2))),
        'maximum_atom_vector_difference_eV_A': float(vectors.max()),
        'marked_pair_projection_delta_eV_A': float(np.dot(df[j] - df[i], direction)),
        'budgets': {'absolute_energy_meV_atom': 2.0, 'force_vector_rms_eV_A': 0.01, 'absolute_pair_projection_eV_A': 0.02},
    }
    result['pass'] = {
        'native_energy': abs(result['native_energy_delta_meV_atom']) <= 2.0,
        'free_energy': abs(result['free_energy_delta_meV_atom']) <= 2.0,
        'force_vector_rms': result['force_vector_rms_eV_A'] <= 0.01,
        'marked_pair_projection': abs(result['marked_pair_projection_delta_eV_A']) <= 0.02,
    }
    return result

archives = {key: verify_archive(key) for key in ARCHIVES}
mesh = pair_metrics(archives['AgSi_k6_sigma010'], archives['AgSi_k8_sigma010'])
dipole = pair_metrics(archives['dipole_off_16A'], archives['dipole_on_16A'])
vacuum = vacuum_pair_metrics(archives['AgSi_k6_sigma010'], archives['AgSi_vacuum26_k6_sigma010'])
mesh_method_match = {k: v for k, v in archives['AgSi_k6_sigma010']['summary']['method'].items() if k != 'kpts'} == {k: v for k, v in archives['AgSi_k8_sigma010']['summary']['method'].items() if k != 'kpts'}
dipole_method_match = {k: v for k, v in archives['dipole_off_16A']['summary']['method'].items() if k != 'poissonsolver'} == {k: v for k, v in archives['dipole_on_16A']['summary']['method'].items() if k != 'poissonsolver'}
if not mesh_method_match or not dipole_method_match:
    raise SystemExit('Matched mesh or dipole method check failed')
mesh['same_method_except_kpts'] = mesh_method_match
dipole['same_method_except_poissonsolver'] = dipole_method_match

proposal_dir = ROOT / 'research/mace-v12-transfer/coordination/reports/window-b/v19_gap_scan_dipole_vacuum_control_proposals_window_b'
proposal_checks = {}
for line in (proposal_dir / 'SHA256SUMS.txt').read_text().splitlines():
    digest, rel = line.split('  ', 1)
    proposal_checks[rel] = sha(proposal_dir / rel) == digest
if not all(proposal_checks.values()):
    raise SystemExit('Vacuum proposal SHA inventory failed')
proposal_manifest = json.loads((proposal_dir / 'proposal_manifest.json').read_text())
parent_record = proposal_manifest['source_record']
gap_root = BASE / 'pbe_interface_v19_gap_scan'
gap_manifest = json.loads((gap_root / 'input_manifest.json').read_text())
gap_parent_record = next(r for r in gap_manifest['records'] if r['label'] == parent_record['label'])
proposal_rows = {}
for record in proposal_manifest['proposals']:
    atoms = read(proposal_dir / record['input'], format='extxyz')
    proposal_rows[record['proposal_id']] = {
        'sha256': sha(proposal_dir / record['input']),
        'cell_z_A': float(atoms.cell.lengths()[2]),
        'coordinate_span_z_A': float(np.ptp(atoms.positions[:, 2])),
        'geometric_empty_length_A': float(atoms.cell.lengths()[2] - np.ptp(atoms.positions[:, 2])),
        'pbc': atoms.pbc.tolist(),
        'launch_enabled': record['launch_enabled'],
        'poissonsolver': record['method']['poissonsolver'],
        'input_sha256': record['input_sha256'],
    }
parent_archive_dir = gap_root / 'calculations' / gap_parent_record['label']
result = {
    'status': 'SLAB_REFERENCE_READINESS_REVIEW',
    'archive_controls': {k: {kk: vv for kk, vv in v.items() if kk not in ('summary', 'atoms')} for k, v in archives.items()},
    'existing_paired_comparisons': {
        'kmesh_6x6x1_to_8x8x1_at_sigma010_16A_parent_method': mesh,
        'vacuum_16A_to_26A_same_k6x6x1_sigma010_all_periodic': vacuum,
        'dipole_off_to_on_at_16A': dipole,
    },
    'vacuum_proposals': {
        'proposal_sha_checks': proposal_checks,
        'source_parent_record': {k: parent_record.get(k) for k in ('label', 'input_sha256', 'pbc', 'poissonsolver', 'method', 'launch_enabled', 'owner_task')},
        'gap_scan_manifest': {'launch_enabled': gap_manifest['launch_enabled'], 'status': gap_manifest.get('status'), 'parent_launch_status': gap_parent_record.get('launch_status'), 'parent_input_sha256': gap_parent_record['input_sha256'], 'parent_method': gap_parent_record['method']},
        'geometry': proposal_rows,
        '26A_parent_archive_exists': (parent_archive_dir / 'verification.json').is_file(),
        '36A_proposal_archive_exists': False,
    },
    'interpretation': [
        'The archive called originalvac has only 16.000267 A geometric empty length (cell z 34.205877 A minus coordinate span 18.205609 A); its dipole ON/OFF pair does not test the 26 A gap-scan vacuum.',
        'The disabled gap-scan parent is a 26.000267 A geometric-vacuum, dipole-ON input. The 26 A OFF proposal changes only the dipole setting and has no DFT label.',
        'The 36 A proposal keeps dipole ON and changes cell z by 10 A while translating all atoms +5 A; no DFT result exists for it.',
        'The sigma=0.1 k6x6x1 to k8x8x1 mesh pair passes the current energy, vector, and marked-pair budgets on the 16 A parent mesh-pilot method. That pair uses its recorded PBC/Poisson settings and is not a matched validation of the 26 A dipole-ON gap-scan method.',
        'The existing 16 A dipole ON/OFF pair has small deltas, but it does not establish vacuum-size convergence or model-vs-DFT interface accuracy.',
        'The 16 A to 26 A all-periodic vacuum pair independently passes the current energy and force budgets. Because it has no dipole correction and uses 3D periodic boundary conditions, it is supporting vacuum evidence rather than the exact 26 A dipole-ON target method.',
        'The mesh and smearing ladders also use fully periodic z. They are not convergence evidence for the nonperiodic-z slab target. The target slab still needs matched 4x4x1/6x6x1/8x8x1 mesh and sigma 0.05/0.10/0.20 controls, followed by a matched vacuum-size comparison before any +/-0.20 A gap points are enabled.',
        'The interface force and separation-force error thresholds remain unassessed; bulk Ti3SiC2 mesh convergence cannot answer them.',
    ],
    'recommended_next_step': 'Run the registered AgSi-owned 26 A dipole-ON parent as a provisional target-boundary reference. Before enabling any +/-0.20 A gap point, perform the nonperiodic-z target-slab 4x4x1/6x6x1/8x8x1 mesh ladder and sigma 0.05/0.10/0.20 controls at fixed geometry, then isolate vacuum size with matched 26 A/36 A calculations under one dipole setting. The existing 26 A dipole-OFF and 36 A dipole-ON proposals change two variables and are not a valid vacuum contrast.',
}
(REPORT / 'audit.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'mesh_pair': mesh, 'vacuum_pair': vacuum, 'dipole_pair': dipole, 'vacuum': {k:v['geometric_empty_length_A'] for k,v in proposal_rows.items()}, '26A_archive': result['vacuum_proposals']['26A_parent_archive_exists'], '36A_archive': result['vacuum_proposals']['36A_proposal_archive_exists']}, indent=2))
