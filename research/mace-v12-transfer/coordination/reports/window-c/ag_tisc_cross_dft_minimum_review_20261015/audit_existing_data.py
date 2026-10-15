#!/usr/bin/env python3
"""Read-only audit of explicitly public reference data; never launches a calculator."""
from pathlib import Path
import gzip
import hashlib
import json
import re
import shutil
import subprocess
from collections import Counter
from datetime import datetime, timezone

import numpy as np
from ase.io import read

OUT = Path(__file__).resolve().parent
REPO = OUT.parents[5]
BUNDLE = REPO / 'research/mace-v12-transfer/exports/cross_potential_dft_reference_bundle_20261009'
E_REPORT = ('research/mace-v12-transfer/coordination/reports/window-e/'
            'low_mismatch_217_dft_lite_20261013/minimal_135atom_model/'
            'potential_development_decision_20261014')

def sha(data):
    return hashlib.sha256(data).hexdigest()

def git(*args):
    return subprocess.check_output(['git', '-C', str(REPO), *args])

def checksum_map(text):
    return dict((line.split(maxsplit=1)[1].lstrip('*'), line.split(maxsplit=1)[0])
                for line in text.splitlines() if line.strip())

def cgroup(name):
    p = Path('/sys/fs/cgroup') / name
    return p.read_text().strip() if p.exists() else None

def main():
    ecommit = git('rev-parse', 'origin/cloud-e-b1-dft').decode().strip()
    paths = git('ls-tree', '-r', '--name-only', ecommit).decode().splitlines()
    handoffs = [p for p in paths if 'AG_TISC_C_DFT_HANDOFF_20261015' in p or p.endswith('/DFT_REQUEST_TO_C.md')]
    esums = checksum_map(git('show', f'{ecommit}:{E_REPORT}/SHA256SUMS.txt').decode())
    ereports = []
    for name, expected in esums.items():
        actual = sha(git('show', f'{ecommit}:{E_REPORT}/{name}'))
        ereports.append({'path': f'{E_REPORT}/{name}', 'expected_sha256': expected,
                         'actual_sha256': actual, 'pass': actual == expected})

    sums = checksum_map((BUNDLE / 'SHA256SUMS.txt').read_text())
    verified_files = {}
    def verify(relative):
        p = BUNDLE / relative
        actual = sha(p.read_bytes())
        expected = sums.get(relative)
        verified_files[relative] = {'actual_sha256': actual, 'expected_sha256': expected,
                                    'pass': expected == actual}
        return expected == actual
    verify('bundle_manifest.json')
    verify('pbe_interface_v19_residual_acquisition/input_manifest.json')
    manifest = json.loads((BUNDLE / 'bundle_manifest.json').read_text())
    im = json.loads((BUNDLE / 'pbe_interface_v19_residual_acquisition/input_manifest.json').read_text())
    parents = {x['label']: x for x in im['records']}
    rows = []
    setup_root = Path('/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups')
    for e in manifest['entries']:
        if e['group'] != 'pbe_interface_v19_residual_acquisition':
            continue
        for key in ('input', 'result_extxyz', 'summary'):
            verify(e[key])
        calcdir = Path(e['summary']).parent
        for filename in ('verification.json', 'gpaw.log', 'SHA256SUMS.txt'):
            verify(str(calcdir / filename))
        summary = json.loads((BUNDLE / e['summary']).read_text())
        verification = json.loads((BUNDLE / calcdir / 'verification.json').read_text())
        log = (BUNDLE / calcdir / 'gpaw.log').read_text()
        artifact_checks = {n: sha((BUNDLE / calcdir / n).read_bytes()) == h
                           for n, h in verification['sha256'].items()}
        source = read(BUNDLE / e['input'])
        result = read(BUNDLE / e['result_extxyz'])
        forces = result.arrays['PW_PBE_forces']
        identities = {
            'input_sha': sha((BUNDLE / e['input']).read_bytes()) == e['input_sha256'] == summary['source_sha256'],
            'symbols': source.get_chemical_symbols() == result.get_chemical_symbols(),
            'positions': bool(np.allclose(source.positions, result.positions, atol=1e-8, rtol=0)),
            'cell': bool(np.allclose(source.cell.array, result.cell.array, atol=1e-8, rtol=0)),
            'pbc': bool(np.array_equal(source.pbc, result.pbc)),
            'atom_ids': bool(np.array_equal(source.arrays['lammps_id'], result.arrays['lammps_id'])),
            'forces_complete_finite': bool(forces.shape == (len(result), 3) and np.isfinite(forces).all()),
            'force_max_matches': bool(np.isclose(np.linalg.norm(forces, axis=1).max(), summary['fmax_eV_A'], atol=1e-7)),
            'native_energy_matches': bool(np.isclose(result.info['PW_PBE_energy_eV'], summary['energy_eV_cell'], atol=1e-10, rtol=0)),
            'scf_summary': summary['scf_converged'] is True,
            'scf_log': bool(re.search(r'Converged in ' + str(summary['scf_iterations']) + r' steps', log)),
            'archive_verifier': verification['status'] == 'PASS' and all(verification['checks'].values()),
            'archive_artifact_sha': all(artifact_checks.values()),
        }
        setup_ids = dict(re.findall(r'^  (Ag|Ti|Si|C):\n.*?^    id:\s+(\w+)', log, re.M | re.S))
        setup_bindings = {}
        for symbol, old_id in setup_ids.items():
            p = setup_root / f'{symbol}.PBE.gz'
            raw = p.read_bytes()
            setup_bindings[symbol] = {'logged_setup_id': old_id,
                'current_setup_xml_md5': hashlib.md5(gzip.decompress(raw)).hexdigest(),
                'current_file_sha256': sha(raw)}
        identities['paw_xml_fingerprints'] = bool(setup_bindings) and all(
            x['logged_setup_id'] == x['current_setup_xml_md5'] for x in setup_bindings.values())
        rows.append({'label': e['label'], 'input_path': e['input'], 'result_path': e['result_extxyz'],
            'natoms': len(result), 'composition': dict(Counter(result.get_chemical_symbols())),
            'cell_A': result.cell.array.tolist(), 'pbc': result.pbc.tolist(),
            'parent_geometry_sha256': parents[e['label']]['parent_geometry_sha256'],
            'original_role': summary['dataset_role'], 'method': summary['method'],
            'native_energy_eV': summary['energy_eV_cell'], 'free_energy_eV': summary['free_energy_eV_cell'],
            'scf_iterations': summary['scf_iterations'], 'elapsed_s': summary['elapsed_s'],
            'force_max_eV_per_A': summary['fmax_eV_A'], 'checks': identities,
            'setup_bindings': setup_bindings, 'result_integrity_pass': all(identities.values()),
            'kpoint_numerical_convergence': 'NOT_ESTABLISHED_BY_THIS_AUDIT',
            'matched_component_energy_ready': False})

    disk = shutil.disk_usage('/workspace')
    memory_max = cgroup('memory.max')
    memory_current = cgroup('memory.current')
    evidence = {'audit_time_UTC': datetime.now(timezone.utc).isoformat(),
        'audit_code_sha256': sha(Path(__file__).read_bytes()),
        'local_base_commit': git('rev-parse', 'HEAD').decode().strip(),
        'e_remote_commit': ecommit, 'requested_handoff_paths_found': handoffs,
        'e_analysis_report_checksums': ereports,
        'public_six_records': rows, 'bundle_files_checked': verified_files,
        'six_records_pass': len(rows) == 6 and all(x['result_integrity_pass'] for x in rows)
                            and all(x['pass'] for x in verified_files.values()),
        'parent_family_count': len({x['parent_geometry_sha256'] for x in rows}),
        'historical_walltime_hours_sum': sum(x['elapsed_s'] for x in rows) / 3600,
        'resources': {'cpu_max': cgroup('cpu.max'), 'memory_max_bytes': memory_max,
            'memory_current_bytes': memory_current, 'memory_headroom_bytes':
                int(memory_max) - int(memory_current) if memory_max and memory_max != 'max' else None,
            'memory_peak_bytes': cgroup('memory.peak'), 'memory_events': cgroup('memory.events'),
            'disk_free_bytes': disk.free, 'pw_x_in_PATH': shutil.which('pw.x'),
            'mpiexec_in_PATH': shutil.which('mpiexec'), 'SCF_started': False},
        'requested_135_matched_triplet_verified': False,
        'new_DFT_execution_authorized': False}
    (OUT / 'AUDIT_EVIDENCE.json').write_text(json.dumps(evidence, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'six_records_pass': evidence['six_records_pass'],
        'parent_families': evidence['parent_family_count'], 'ecommit': ecommit,
        'requested_handoff_found': bool(handoffs), 'resources': evidence['resources'],
        'records': [{'label': r['label'], 'atoms':r['natoms'], 'pass':r['result_integrity_pass'],
                     'checks':r['checks'], 'hours':r['elapsed_s']/3600} for r in rows]}, indent=2))

if __name__ == '__main__':
    main()
