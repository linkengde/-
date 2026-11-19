#!/workspace/venvs/dft-labels/bin/python
"""Prepare isolated v12m refit from v12l with a DFT-labelled Ag-C-C collapse precursor.

The established 128/12 train/validation split, 21 existing local force frames, and Si01 center holdout are inherited unchanged from v12l. Local crop edge forces are anchored to POD v10;
only the two central Ti-C force rows in each crop use DFT targets.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write

HERE = Path(__file__).resolve().parent
OUT = HERE / 'v12m_agctriad_candidate'
V12L = HERE / 'v12l_ticpair_candidate'
ACTIVE = HERE.parent / 'target_200x200x300_qsgs_material200' / 'active_tic_hotspot_expansion'
V10 = HERE.parent / 'final_v10' / 'POD'
LMP = Path('/workspace/software/lammps-pod-snap-clean/build/lmp')
ELEMENTS = ['Ag', 'C', 'Si', 'Ti']
MASS = {'Ag': 107.8682, 'C': 12.011, 'Si': 28.085, 'Ti': 47.867}
LABELS = [
    ('agc_triad_v12l_t28fs_step18', 'agc_triad_v12l_t28fs_step18_pw350_gamma', [414900, 198354, 203868]),
]

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def run_v10_pod(atoms: Atoms, work: Path) -> np.ndarray:
    work.mkdir(parents=True)
    data = work / 'local.data'
    write(data, atoms, format='lammps-data', atom_style='atomic', units='metal', specorder=ELEMENTS)
    masses = '\n'.join(f'mass {i} {MASS[s]}' for i, s in enumerate(ELEMENTS, 1))
    inp = work / 'forces.in'
    inp.write_text(f'''clear
units metal
atom_style atomic
boundary p p p
read_data {data}
{masses}
pair_style pod
pair_coeff * * {V10/'AgTiSiC_POD.pod'} {V10/'pod_coefficients.pod'} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
run 0
write_dump all custom forces.dump id type x y z fx fy fz modify sort id
''')
    subprocess.run([str(LMP), '-in', str(inp), '-log', str(work/'lammps.log'), '-screen', 'none'],
                   cwd=work, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    lines = (work/'forces.dump').read_text().splitlines()
    h = next(i for i, line in enumerate(lines) if line.startswith('ITEM: ATOMS '))
    cols = lines[h].split()[2:]
    rows = np.loadtxt(lines[h+1:], usecols=[cols.index(k) for k in ('id', 'fx', 'fy', 'fz')],
                      dtype=float, ndmin=2)
    rows = rows[np.argsort(rows[:, 0])]
    if len(rows) != len(atoms) or not np.isfinite(rows[:, 1:]).all():
        raise RuntimeError('Invalid POD v10 baseline force output')
    return rows[:, 1:4]


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f'Refusing to overwrite {OUT}')
    base_local = V12L / 'dataset/local_train/local_dft_targets.xyz'
    label_paths = [ACTIVE / f'labels/{stem}.extxyz' for _, stem, _ in LABELS]
    for path in (base_local, V10/'AgTiSiC_POD.pod', V10/'pod_coefficients.pod', LMP,
                 V12L/'dataset/validation', V12L/'POD/fitpod.data', V12L/'SNAP/fitsnap.in',
                 *label_paths):
        if not path.exists():
            raise FileNotFoundError(path)

    sources = read(base_local, index=':')
    if len(sources) != 21:
        raise RuntimeError(f'Expected 21 inherited local frames, found {len(sources)}')
    # Copy unchanged baseline training/validation roles. Add one targeted Ag-C-C crop.
    shutil.copytree(V12L/'dataset/train', OUT/'dataset/train',
                    ignore=shutil.ignore_patterns('local_dft_targets.xyz'))
    shutil.copytree(V12L/'dataset/validation', OUT/'dataset/validation')
    shutil.copytree(V12L/'SNAP/data', OUT/'SNAP/data',
                    ignore=shutil.ignore_patterns('Local.xyz'))
    (OUT/'dataset/local_train').mkdir(parents=True)
    (OUT/'POD').mkdir(parents=True)
    for name in ('AgTiSiC_POD_v12l.pod', 'fitpod.data', 'fitpod.in'):
        shutil.copy2(V12L/'POD'/name, OUT/'POD'/name)
    for name in ('fitsnap.in',):
        shutil.copy2(V12L/'SNAP'/name, OUT/'SNAP'/name)
    shutil.copy2(V12L/'holdout_si_hotspot_104359.extxyz',
                 OUT/'holdout_si_hotspot_104359.extxyz')
    holdout_data = HERE/'holdout_label_evaluation/POD_v12f_si01/local.data'
    holdout_data_dst = OUT/'holdout_label_evaluation/POD_v12f_si01/local.data'
    holdout_data_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(holdout_data, holdout_data_dst)

    added = []
    fit_frames = []
    for (state, stem, expected_centers), path in zip(LABELS, label_paths):
        log = ACTIVE / f'labels/{stem}.txt'
        record = json.loads((ACTIVE / f'labels/{stem}.json').read_text())
        if record.get('status') != 'CONVERGED' or 'Converged in ' not in log.read_text(errors='replace'):
            raise RuntimeError(f'{state} label is not SCF-converged')
        if record.get('center_atom_ids') != expected_centers or record.get('DFT_target_rows') != 3:
            raise RuntimeError(f'Unexpected Ag-C-C target IDs or target row count for {state}')
        dft = read(path)
        ids = np.asarray(dft.arrays['source_atom_id'], dtype=int)
        target = np.asarray(dft.arrays['is_DFT_force_target'], dtype=bool)
        dft_forces = np.asarray(dft.arrays['REF_forces'], dtype=float)
        if int(target.sum()) != 3 or not np.isfinite(dft_forces).all():
            raise RuntimeError(f'Invalid DFT force data for {state}')
        if sorted(ids[target].tolist()) != sorted(expected_centers):
            raise RuntimeError(f'Target atom mapping mismatch for {state}')
        baseline = run_v10_pod(dft, OUT/'local_pod_v10_forces'/state)
        fit_forces = baseline.copy()
        fit_forces[target] = dft_forces[target]
        rich = Atoms(symbols=dft.get_chemical_symbols(), positions=dft.positions,
                     cell=dft.cell, pbc=dft.pbc)
        rich.new_array('forces', fit_forces)
        rich.new_array('source_atom_id', ids.astype(np.int32))
        rich.new_array('DFT_forces', dft_forces)
        rich.new_array('POD_v10_forces', baseline)
        rich.new_array('is_DFT_force_target', target.astype(np.int8))
        rich.info.update({'energy': float(dft.info['REF_energy']),
                          'config_type': f'POD_v12m_full200_AgCC_{state}',
                          'center_atom_ids': json.dumps(expected_centers),
                          'DFT_method': str(dft.info['DFT_method']),
                          'DFT_label_sha256': sha256(path),
                          'role': 'training_candidate; SI01 center 104359 remains held out'})
        sources.append(rich)
        clean = Atoms(symbols=rich.get_chemical_symbols(), positions=rich.positions,
                      cell=rich.cell, pbc=rich.pbc)
        clean.new_array('forces', fit_forces.copy())
        clean.info['energy'] = float(dft.info['REF_energy'])
        clean.info['config_type'] = rich.info['config_type']
        fit_frames.append(clean)
        added.append({'state': state, 'natoms': len(dft), 'center_pair_distances_A': record['center_pair_distances_A'],
                      'target_atom_ids': expected_centers,
                      'DFT_target_force_norms_eV_A': np.linalg.norm(dft_forces[target], axis=1).tolist(),
                      'source_model_center_force_errors_eV_A': record.get('center_POD_v12l_force_errors_eV_A'),
                      'DFT_label_sha256': sha256(path), 'SCF_iterations': record['iterations']})

    rich_path = OUT/'dataset/local_train/local_dft_targets.xyz'
    write(rich_path, sources, format='extxyz')
    # FitPOD's reader uses positional columns, so fitter-only extxyz files
    # preserve strict species,pos,forces ordering.
    fit_frames = []
    for src in sources:
        if 'forces' in src.arrays:
            force = np.asarray(src.arrays['forces'], dtype=float)
        else:
            force = np.asarray(src.get_forces(apply_constraint=False), dtype=float)
        clean = Atoms(symbols=src.get_chemical_symbols(), positions=src.positions,
                      cell=src.cell, pbc=src.pbc)
        clean.new_array('forces', force.copy())
        if 'energy' in src.info:
            clean.info['energy'] = float(src.info['energy'])
        else:
            clean.info['energy'] = float(src.get_potential_energy())
        clean.info['config_type'] = str(src.info['config_type'])
        fit_frames.append(clean)
    write(OUT/'dataset/train/local_dft_targets.xyz', fit_frames, format='extxyz')
    write(OUT/'SNAP/data/Local.xyz', fit_frames, format='extxyz')
    audit_files = [OUT/'dataset/train/local_dft_targets.xyz', OUT/'SNAP/data/Local.xyz']
    for path in audit_files:
        check = read(path, index=':')
        if len(check) != len(fit_frames):
            raise RuntimeError(f'Frame count mismatch in {path}')
        for i, (expected, actual) in enumerate(zip(fit_frames, check)):
            if not np.array_equal(expected.numbers, actual.numbers):
                raise RuntimeError(f'Species order changed in {path}, frame {i}')
            if not np.allclose(expected.positions, actual.positions, atol=0, rtol=0):
                raise RuntimeError(f'Positions changed in {path}, frame {i}')
            if not np.allclose(expected.arrays['forces'], actual.get_forces(), atol=1e-12, rtol=0):
                raise RuntimeError(f'Forces changed in {path}, frame {i}')
            if float(expected.info['energy']) != float(actual.get_potential_energy()):
                raise RuntimeError(f'Energy changed in {path}, frame {i}')
        if 'Properties=species:S:1:pos:R:3:forces:R:3' not in path.read_text().splitlines()[1]:
            raise RuntimeError(f'Unexpected force-column order in {path}')
    (OUT/'local_fit_input_audit.json').write_text(json.dumps({
        'status': 'PASS', 'local_frames': len(fit_frames),
        'fitter_files': [str(p) for p in audit_files],
        'property_order': 'species, position[3], forces[3]',
        'new_DFT_target_rows': 3,
        'total_DFT_target_force_rows_in_local_frames': 17,
        'note': 'For the new Ag-C-C crop, only the three target atom forces are DFT targets; all other crop forces are POD v10 anchors.'
    }, indent=2, ensure_ascii=False)+'\n')

    desc = OUT/'POD/AgTiSiC_POD_v12m.pod'
    shutil.copy2(V12L/'POD/AgTiSiC_POD_v12l.pod', desc)
    cfg = (OUT/'POD/fitpod.data').read_text()
    cfg = cfg.replace(str(V12L/'dataset/train'), str(OUT/'dataset/train'))
    cfg = cfg.replace(str(V12L/'dataset/validation'), str(OUT/'dataset/validation'))
    cfg = cfg.replace('basename_for_output_files pod_v12l', 'basename_for_output_files pod_v12m')
    (OUT/'POD/fitpod.data').write_text(cfg)
    (OUT/'POD/fitpod.in').write_text('fitpod AgTiSiC_POD_v12m.pod fitpod.data\n')
    snap = (OUT/'SNAP/fitsnap.in').read_text()
    snap = snap.replace(str(V12L/'dataset/train'), str(OUT/'dataset/train'))
    snap = snap.replace(str(V12L/'dataset/validation'), str(OUT/'dataset/validation'))
    snap = snap.replace('SNAP_v12l_metrics.md', 'SNAP_v12m_metrics.md')
    snap = snap.replace('AgTiSiC_SNAP_v12l', 'AgTiSiC_SNAP_v12m')
    snap_lines = [line for line in snap.splitlines() if not line.strip().startswith('bnormflag')]
    chem_index = next((i for i, line in enumerate(snap_lines) if line.strip() == 'chemflag = 1'), None)
    if chem_index is None:
        raise RuntimeError('SNAP config is missing chemflag = 1')
    snap_lines.insert(chem_index + 1, 'bnormflag = 1')
    snap = '\n'.join(snap_lines) + '\n'
    if snap.count('bnormflag = 1') != 1:
        raise RuntimeError('SNAP config must contain exactly one bnormflag = 1')
    (OUT/'SNAP/fitsnap.in').write_text(snap)

    manifest = json.loads((V12L/'local_training_manifest.json').read_text())
    manifest.update({'status': 'prepared; waiting for independent fit and validation',
                     'candidate': 'v12m_full200_Ag-C-C_cold-stability_refit',
                     'local_force_frames': len(sources),
                     'si01_holdout_center_atom_id': 104359,
                     'holdout_excluded_from_training': True,
                     'new_tic_hotspot_labels': list(manifest.get('new_tic_hotspot_labels', [])) + added,
                     'fit_input_property_order': 'species:S:1,pos:R:3,forces:R:3',
                     'SNAP_flags': {'chemflag': 1, 'bnormflag': 1},
                     'total_DFT_target_force_rows': 17,
                     'local_force_semantics': 'DFT forces on Ag-C-C target atoms; POD v10 force anchors on crop edge atoms; local force group weight=10.0'})
    manifest['source_labels'] = list(dict.fromkeys(
        list(manifest.get('source_labels', [])) + [str(p) for p in label_paths]))
    manifest['total_fitpod_training_configurations'] = 128 + len(sources)
    (OUT/'local_training_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+'\n')
    (OUT/'run_v12m_fits.sh').write_text('''#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
(cd "$root/POD" && "/workspace/software/lammps-pod-snap-clean/build/lmp" -in fitpod.in -log fitpod_v12m.log)
test -s "$root/POD/pod_v12m_coefficients.pod"
test -s "$root/POD/pod_v12m_test_errors.pod"
(cd "$root/SNAP" && "/workspace/venvs/fitsnap/bin/python" /workspace/potential-fit/round7_pod_snap/run_fitsnap_py312.py --overwrite fitsnap.in > fitsnap_v12m.log 2>&1)
test -s "$root/SNAP/AgTiSiC_SNAP_v12m.snapcoeff"
test -s "$root/SNAP/AgTiSiC_SNAP_v12m.snapparam"
test -s "$root/SNAP/SNAP_v12m_metrics.md"
''')
    (OUT/'run_v12m_fits.sh').chmod(0o755)
    evaluator = (V12L/'evaluate_v12l_si01.py').read_text()
    evaluator = evaluator.replace('v12l', 'v12m')
    (OUT/'evaluate_v12m_si01.py').write_text(evaluator)
    (OUT/'evaluate_v12m_si01.py').chmod(0o755)
    full_root = HERE.parent/'target_200x200x300_qsgs_material200'
    cold_full = (full_root/'run_cold_v12l_initial_458k.py').read_text()
    cold_full = cold_full.replace('v12l_ticpair_candidate', 'v12m_agctriad_candidate').replace('v12l', 'v12m')
    (full_root/'run_cold_v12m_initial_458k.py').write_text(cold_full)
    (full_root/'run_cold_v12m_initial_458k.py').chmod(0o755)
    static_full = (full_root/'run_static_v12l_458k.py').read_text()
    static_full = static_full.replace('v12l_ticpair_candidate', 'v12m_agctriad_candidate').replace('v12l', 'v12m')
    (full_root/'run_static_v12m_458k.py').write_text(static_full)
    (full_root/'run_static_v12m_458k.py').chmod(0o755)
    static_snap = (full_root/'run_static_snap_v12l_458k.py').read_text()
    static_snap = static_snap.replace('v12l_ticpair_candidate', 'v12m_agctriad_candidate').replace('v12l', 'v12m')
    (full_root/'run_static_snap_v12m_458k.py').write_text(static_snap)
    (full_root/'run_static_snap_v12m_458k.py').chmod(0o755)
    print(json.dumps({'status': 'PREPARED', 'output': str(OUT), 'baseline_train': 128,
                      'local_train_frames': len(sources), 'validation': 12,
                      'SI01_holdout_center_atom_id': 104359, 'new_labels': added},
                     indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
