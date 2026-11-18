#!/workspace/venvs/dft-labels/bin/python
"""Prepare an isolated v12g refit with Si hotspot 04; preserve v12f and SI01."""
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
OUT = HERE / 'v12g_si04_candidate'
ACTIVE = HERE / 'active_si_hotspot_expansion'
LABEL = ACTIVE / 'labels/si_hotspot_04_pw350_gamma.extxyz'
LABEL_LOG = ACTIVE / 'labels/si_hotspot_04_pw350_gamma.txt'
V12F = HERE
V10 = HERE.parent / 'final_v10/POD'
LMP = Path('/workspace/software/lammps-pod-snap-clean/build/lmp')
ELEMENTS = ['Ag', 'C', 'Si', 'Ti']
MASS = {'Ag': 107.8682, 'C': 12.011, 'Si': 28.085, 'Ti': 47.867}
CENTER_ID = 244669


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
    use = [cols.index(k) for k in ('id', 'fx', 'fy', 'fz')]
    rows = np.loadtxt(lines[h+1:], usecols=use, dtype=float, ndmin=2)
    rows = rows[np.argsort(rows[:, 0])]
    if len(rows) != len(atoms) or not np.isfinite(rows[:, 1:]).all():
        raise RuntimeError('Invalid v10 baseline forces for Si hotspot 04')
    return rows[:, 1:4]


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f'Refusing to overwrite {OUT}')
    for p in (LABEL, LABEL_LOG, V10/'AgTiSiC_POD.pod', V10/'pod_coefficients.pod', LMP):
        if not p.is_file():
            raise FileNotFoundError(p)
    log = LABEL_LOG.read_text(errors='replace')
    if 'Converged in ' not in log:
        raise RuntimeError('GPAW log does not confirm SCF convergence; candidate remains excluded')
    dft = read(LABEL)
    if int(dft.info.get('center_atom_id', -1)) != CENTER_ID:
        raise RuntimeError('Unexpected candidate center ID')
    if dft.info.get('DFT_method') != 'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV':
        raise RuntimeError('Unexpected DFT method or cutoff')
    ids = np.asarray(dft.arrays['source_atom_id'], dtype=int)
    match = np.flatnonzero(ids == CENTER_ID)
    forces = np.asarray(dft.arrays['REF_forces'], dtype=float)
    if len(match) != 1 or forces.shape != (len(dft), 3) or not np.isfinite(forces).all():
        raise RuntimeError('Invalid center ID mapping or non-finite DFT forces')
    center = int(match[0])
    base_forces = run_v10_pod(dft, OUT/'local_pod_v10_forces/si_hotspot_04')

    shutil.copytree(V12F/'dataset/train', OUT/'dataset/train',
                    ignore=shutil.ignore_patterns('local_dft_targets.xyz'))
    shutil.copytree(V12F/'dataset/validation', OUT/'dataset/validation')
    (OUT/'dataset/local_train').mkdir(parents=True)
    shutil.copytree(V12F/'SNAP/data', OUT/'SNAP/data',
                    ignore=shutil.ignore_patterns('Local.xyz'))
    shutil.copytree(V12F/'POD', OUT/'POD',
                    ignore=shutil.ignore_patterns('pod_v12f_*', 'stdout_v12f.log'))
    shutil.copytree(V12F/'blind_round7', OUT/'blind_round7')
    source_label = OUT/'source_label_si04'
    source_label.mkdir(parents=True)
    for name in ('si_hotspot_04_candidate.extxyz', 'si_hotspot_04_candidate.json'):
        shutil.copy2(ACTIVE/name, source_label/name)
    for name in ('si_hotspot_04_pw350_gamma.extxyz', 'si_hotspot_04_pw350_gamma.json',
                 'si_hotspot_04_pw350_gamma.txt'):
        shutil.copy2(ACTIVE/'labels'/name, source_label/name)
    for name in ('holdout_si_hotspot_104359.extxyz', 'local_training_manifest.json',
                 'holdout_label_evaluation/summary.json'):
        dest = OUT/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(V12F/name, dest)

    baseline_path = V12F/'dataset/local_train/local_dft_targets.xyz'
    frames = read(baseline_path, index=':')
    fit_forces = base_forces.copy()
    fit_forces[center] = forces[center]
    fit = Atoms(symbols=dft.get_chemical_symbols(), positions=dft.positions,
                cell=dft.cell, pbc=dft.pbc)
    fit.new_array('forces', fit_forces)
    fit.new_array('source_atom_id', ids.astype(np.int32))
    fit.new_array('DFT_forces', forces)
    fit.new_array('POD_v10_forces', base_forces)
    mask = np.zeros(len(dft), dtype=np.int8)
    mask[center] = 1
    fit.new_array('is_DFT_force_target', mask)
    fit.info['energy'] = float(dft.info['REF_energy'])
    fit.info['config_type'] = 'POD_v12g_Si_hotspot_04'
    frames.append(fit)
    local = OUT/'dataset/local_train/local_dft_targets.xyz'
    write(local, frames, format='extxyz')
    # FitPOD's extxyz reader consumes the first six per-atom numeric fields as
    # position[3] + force[3] and does not map custom Properties by name. Keep
    # the rich source/audit arrays above, but write minimal fitter copies with
    # forces immediately following positions so IDs/masks cannot be read as Fx.
    fitter_frames = []
    for source in frames:
        source_forces = np.asarray(source.arrays['forces'] if 'forces' in source.arrays
                                   else source.get_forces(apply_constraint=False), dtype=float)
        if source_forces.shape != (len(source), 3) or not np.isfinite(source_forces).all():
            raise RuntimeError('Invalid force rows in local training frame')
        clean = Atoms(symbols=source.get_chemical_symbols(), positions=source.positions,
                      cell=source.cell, pbc=source.pbc)
        clean.new_array('forces', source_forces.copy())
        energy = source.info.get('energy')
        clean.info['energy'] = float(energy if energy is not None else source.get_potential_energy())
        clean.info['config_type'] = str(source.info.get('config_type', 'local_force_target'))
        fitter_frames.append(clean)
    write(OUT/'dataset/train/local_dft_targets.xyz', fitter_frames, format='extxyz')
    write(OUT/'SNAP/data/Local.xyz', fitter_frames, format='extxyz')
    shutil.copy2(HERE/'sanitize_v12g_fit_data.py', OUT/'sanitize_local_fit_data.py')

    pod_desc = OUT/'POD/AgTiSiC_POD_v12g.pod'
    shutil.copy2(V12F/'POD/AgTiSiC_POD_v12f.pod', pod_desc)
    cfg = (OUT/'POD/fitpod.data').read_text()
    cfg = cfg.replace(str(V12F/'dataset/train'), str(OUT/'dataset/train'))
    cfg = cfg.replace(str(V12F/'dataset/validation'), str(OUT/'dataset/validation'))
    cfg = cfg.replace('basename_for_output_files pod_v12f', 'basename_for_output_files pod_v12g')
    (OUT/'POD/fitpod.data').write_text(cfg)
    (OUT/'POD/fitpod.in').write_text('fitpod AgTiSiC_POD_v12g.pod fitpod.data\n')

    snap = (V12F/'SNAP/fitsnap.in').read_text()
    snap = snap.replace('SNAP_v12f_metrics.md', 'SNAP_v12g_metrics.md')
    snap = snap.replace('AgTiSiC_SNAP_v12f', 'AgTiSiC_SNAP_v12g')
    (OUT/'SNAP/fitsnap.in').write_text(snap)

    manifest = json.loads((V12F/'local_training_manifest.json').read_text())
    manifest['status'] = 'prepared; independent Si hotspot 04 label audited for SCF convergence and finite forces'
    manifest['local_force_frames'] = len(frames)
    manifest['new_hotspot_labels_in_training'] = [70588, 83874, 6223, CENTER_ID]
    manifest['new_hotspot_labels_held_out'] = [104359]
    manifest['source_labels'].append(str(LABEL))
    manifest['new_si_hotspot_04'] = {
        'center_atom_id': CENTER_ID,
        'target_atom_ids': [CENTER_ID],
        'force_target_rows': int(mask.sum()),
        'dft_center_force_eV_A': forces[center].tolist(),
        'dft_force_norm_eV_A': float(np.linalg.norm(forces[center])),
        'DFT_label_sha256': sha256(LABEL),
        'role': 'training_candidate; SI01 remains development holdout'
    }
    (OUT/'local_training_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps({'status':'PREPARED','output':str(OUT),'train_frames':128,
        'local_force_frames':len(frames),'validation_frames':12,'center_atom_id':CENTER_ID,
        'target_rows':int(mask.sum()),'si01_holdout_id':104359,
        'label_sha256':sha256(LABEL)},indent=2))


if __name__ == '__main__':
    main()
