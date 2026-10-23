#!/usr/bin/env python3
from pathlib import Path
import csv, hashlib, json, os, shutil, subprocess
import numpy as np
from ase.io import read

ROOT = Path('/tmp/e-tsc69-fresh/research/mace-v12-transfer/coordination/reports/window-e')
PKG = ROOT / 'tsc_surface_10m11_dft_diagnostic_20261023'
TASK = 'TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_RESP_TIZ020'
INPUT = ROOT / 'tsc_surface_geometry_dev_candidate_10m11_20261023/inputs' / f'{TASK}.extxyz'
RESULT_DIR = PKG / 'results' / TASK
OUT = RESULT_DIR / 'tscbop_run0' / 'output'
POTENTIAL = ROOT / 'ti3sic2_backbone_e5_force_jacobian_20261014/potentials/TSCBOP.source.tersoff'
POT_SHA = '729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7'
LMP = Path('/workspace/.venvs/lammps-run0/bin/lmp')
TYPE = {'Ti': 1, 'Si': 2, 'C': 3}
MASS = {'Ti': 47.867, 'Si': 28.085, 'C': 12.011}

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

def geometry_sha(atoms):
    h = hashlib.sha256()
    h.update(np.asarray(atoms.numbers, dtype='<i4').tobytes())
    h.update(np.asarray(atoms.cell.array, dtype='<f8').tobytes())
    h.update(np.asarray(atoms.positions, dtype='<f8').tobytes())
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays['source_atom_id'], dtype='<i8').tobytes())
    return h.hexdigest()

def read_dump(path):
    lines = path.read_text().splitlines()
    sec = next(i for i, line in enumerate(lines) if line.startswith('ITEM: ATOMS'))
    cols = lines[sec].split()[2:]
    idx = {k:i for i,k in enumerate(cols)}
    type_to_element = {v:k for k,v in TYPE.items()}
    rows = {}
    for line in lines[sec+1:]:
        if not line.strip():
            continue
        r = line.split()
        atom_id = int(r[idx['id']])
        elem = type_to_element[int(r[idx['type']])]
        xyz = np.array([float(r[idx[k]]) for k in ('x','y','z')])
        frc = np.array([float(r[idx[k]]) for k in ('fx','fy','fz')])
        if not np.isfinite(xyz).all() or not np.isfinite(frc).all():
            raise ValueError(f'nonfinite LAMMPS record for source ID {atom_id}')
        rows[atom_id] = (elem, xyz, frc)
    return rows

def metrics(err, ids, mask):
    e = err[mask]
    norms = np.linalg.norm(e, axis=1)
    local = np.flatnonzero(mask)
    return {
        'atom_count': int(mask.sum()),
        'component_rmse_eV_A': float(np.sqrt(np.mean(e**2))),
        'vector_rmse_eV_A': float(np.sqrt(np.mean(np.sum(e*e, axis=1)))),
        'p95_vector_error_eV_A': float(np.percentile(norms, 95)),
        'max_vector_error_eV_A': float(norms.max()),
        'max_error_source_atom_id': int(ids[local[np.argmax(norms)]])
    }

result = json.loads((RESULT_DIR/'RESULT.json').read_text())
status = json.loads((RESULT_DIR/'RUN_STATUS.json').read_text())
assert result['status'] == 'SUCCESS' and result['scf_converged'] is True
assert status['status'] == 'SUCCESS' and result['scf_converged'] is True
assert sha(INPUT) == result['input_sha256']
assert sha(POTENTIAL) == POT_SHA
assert LMP.is_file()
atoms = read(INPUT, format='extxyz')
labelled = read(RESULT_DIR/'labelled.extxyz', format='extxyz')
assert geometry_sha(atoms) == result['geometry_sha256']
assert geometry_sha(labelled) == result['geometry_sha256']
for field in ('numbers', 'positions'):
    assert np.array_equal(getattr(atoms, field), getattr(labelled, field)), field
assert np.array_equal(atoms.cell.array, labelled.cell.array)
assert np.array_equal(atoms.pbc, labelled.pbc)
ids = np.asarray(atoms.arrays['source_atom_id'], dtype=np.int64)
assert ids.shape == (69,) and len(set(ids.tolist())) == 69 and ids.min() > 0
symbols = atoms.get_chemical_symbols()
assert len(atoms) == 69 and {e: symbols.count(e) for e in set(symbols)} == {'Ti':33,'Si':12,'C':24}
dft = np.asarray(result['forces_eV_per_angstrom'], dtype=float)
assert dft.shape == (69, 3) and np.isfinite(dft).all()
OUT.mkdir(parents=True, exist_ok=True)

# Rotate only the coordinate frame into LAMMPS restricted-triclinic form.
# Fractional coordinates, cell metrics, periodicity and atom IDs are preserved.
cell = np.asarray(atoms.cell.array, dtype=float)
a, b, c = cell
ex = a / np.linalg.norm(a)
bperp = b - np.dot(b, ex)*ex
ey = bperp / np.linalg.norm(bperp)
ez = np.cross(ex, ey)
Q = np.column_stack((ex, ey, ez))
cell_lmp = cell @ Q
if np.linalg.det(cell_lmp) <= 0:
    raise ValueError('rotation produced nonpositive cell volume')
if max(abs(cell_lmp[0,1]), abs(cell_lmp[0,2]), abs(cell_lmp[1,2])) > 1e-10:
    raise ValueError(f'failed to construct restricted-triclinic cell: {cell_lmp}')
frac = np.asarray(atoms.positions) @ np.linalg.inv(cell)
frac = frac - np.floor(frac + 1e-12)
xyz = frac @ cell_lmp
lx, ly, lz = cell_lmp[0,0], cell_lmp[1,1], cell_lmp[2,2]
xy, xz, yz = cell_lmp[1,0], cell_lmp[2,0], cell_lmp[2,1]

with (OUT/'atoms.data').open('w') as f:
    f.write('Frozen TSCBOP-only source geometry; coordinate frame rotated, IDs retained\n\n')
    f.write(f'{len(atoms)} atoms\n3 atom types\n\n')
    f.write(f'0 {lx:.16g} xlo xhi\n0 {ly:.16g} ylo yhi\n0 {lz:.16g} zlo zhi\n')
    f.write(f'{xy:.16g} {xz:.16g} {yz:.16g} xy xz yz\n\nMasses\n\n')
    for el in ('Ti','Si','C'):
        f.write(f'{TYPE[el]} {MASS[el]}\n')
    f.write('\nAtoms # atomic\n\n')
    for atom_id, el, r in zip(ids, symbols, xyz):
        f.write(f'{int(atom_id)} {TYPE[el]} {r[0]:.16g} {r[1]:.16g} {r[2]:.16g}\n')
with (OUT/'coordinate_map.csv').open('w', newline='') as f:
    w = csv.writer(f, lineterminator='\n')
    w.writerow(['source_atom_id','element','x_source_A','y_source_A','z_source_A','x_lammps_A','y_lammps_A','z_lammps_A'])
    for atom_id, el, p, q in zip(ids, symbols, atoms.positions, xyz):
        w.writerow([int(atom_id), el, *map(float,p), *map(float,q)])
shutil.copy2(POTENTIAL, OUT/'TSCBOP.source.tersoff')
assert sha(OUT/'TSCBOP.source.tersoff') == POT_SHA
(OUT/'run.in').write_text(f'''clear\nunits metal\natom_style atomic\nboundary p p p\nread_data {(OUT/'atoms.data').resolve()}\npair_style tersoff\npair_coeff * * {(OUT/'TSCBOP.source.tersoff').resolve()} Ti Si C\nneighbor 2.0 bin\nneigh_modify delay 0 every 1 check yes\nthermo 1\nthermo_style custom step atoms pe\ndump forces all custom 1 forces.dump id type x y z fx fy fz\ndump_modify forces sort id format line "%d %d %.16g %.16g %.16g %.16g %.16g %.16g"\nrun 0 post no\nundump forces\n''')
(OUT/'RUN0_MANIFEST.json').write_text(json.dumps({
    'scope':'Frozen-coordinate LAMMPS run 0 only; no minimization or time integration.',
    'structure_id':TASK,'input_sha256':sha(INPUT),'geometry_sha256':result['geometry_sha256'],
    'atom_count':len(atoms),'composition':{'Ti':33,'Si':12,'C':24},'source_ids_retained':True,
    'coordinate_frame_rotation_only':True,'pbc':[bool(x) for x in atoms.pbc],
    'potential_sha256':POT_SHA,'pair_style':'tersoff','pair_coeff_mapping':['Ti','Si','C'],
    'lammps_executable':str(LMP),'units':'metal; force eV/Angstrom','run_steps':0,
    'omp_threads':1,'dft_scf_steps':result['scf_iterations']},indent=2)+'\n')
env = os.environ.copy(); env['OMP_NUM_THREADS']='1'; env['OMP_DYNAMIC']='FALSE'
if not (OUT/'forces.dump').is_file():
    p = subprocess.run([str(LMP), '-in', 'run.in', '-log', 'lammps.log', '-screen', 'screen.log', '-nocite'], cwd=OUT, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (OUT/'launcher.log').write_text(p.stdout)
    if p.returncode != 0 or not (OUT/'forces.dump').is_file():
        raise RuntimeError(f'LAMMPS run0 failed, exit={p.returncode}; see launcher.log/screen.log')
if 'run 0 post no' not in (OUT/'run.in').read_text() or 'for 0 steps' not in (OUT/'screen.log').read_text(errors='replace'):
    raise RuntimeError('LAMMPS output/input does not confirm run 0')
lmp = read_dump(OUT/'forces.dump')
if set(lmp) != set(ids.tolist()):
    raise ValueError('LAMMPS atom IDs do not match source structure')
force_lmp = np.zeros_like(dft)
maxcoord = 0.0
for j, atom_id in enumerate(ids):
    el, coord, frc = lmp[int(atom_id)]
    if el != symbols[j]:
        raise ValueError(f'element mismatch for source ID {atom_id}: {el} != {symbols[j]}')
    delta_frac = (coord - xyz[j]) @ np.linalg.inv(cell_lmp)
    delta_frac -= np.rint(delta_frac)
    delta_cart = delta_frac @ cell_lmp
    maxcoord = max(maxcoord, float(np.linalg.norm(delta_cart)))
    if np.linalg.norm(delta_cart) > 1e-8:
        raise ValueError(f'periodic coordinate mismatch after run 0 at source ID {atom_id}')
    force_lmp[j] = frc
err = force_lmp - dft
rows = []
for j, atom_id in enumerate(ids):
    rows.append({
      'source_atom_id':int(atom_id),'element':symbols[j],
      'x_A':float(atoms.positions[j,0]),'y_A':float(atoms.positions[j,1]),'z_A':float(atoms.positions[j,2]),
      'dft_fx_eV_A':float(dft[j,0]),'dft_fy_eV_A':float(dft[j,1]),'dft_fz_eV_A':float(dft[j,2]),
      'tscbop_fx_eV_A':float(force_lmp[j,0]),'tscbop_fy_eV_A':float(force_lmp[j,1]),'tscbop_fz_eV_A':float(force_lmp[j,2]),
      'error_fx_eV_A':float(err[j,0]),'error_fy_eV_A':float(err[j,1]),'error_fz_eV_A':float(err[j,2]),
      'force_error_norm_eV_A':float(np.linalg.norm(err[j]))})
with (OUT/'FORCE_COMPARISON.csv').open('w', newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
by_element={el:metrics(err,ids,np.asarray(symbols)==el) for el in ('Ti','Si','C')}
summary={
 'structure_id':TASK,'dft_scf_converged':True,'scf_iterations':result['scf_iterations'],
 'dft_energy_eV':result['energy_eV'],'dft_force_consistent_free_energy_eV':result['force_consistent_free_energy_eV'],
 'dft_input_sha256':sha(INPUT),'dft_geometry_sha256':result['geometry_sha256'],
 'tscbop_sha256':POT_SHA,'lammps_version':'22 Jul 2025 Update 4','units':'metal; eV, Angstrom, eV/Angstrom',
 'atom_count':len(atoms),'coordinate_frame_rotation_only':True,'max_coordinate_roundtrip_error_A':maxcoord,
 'all_atoms':metrics(err,ids,np.ones(len(ids),dtype=bool)),'by_element':by_element,
 'interpretation':'Same frozen structure and source atom IDs. The DFT and TSCBOP-only force comparison is a same-parent-family diagnostic, not an independent blind validation. No fitting or role assignment.'}
(OUT/'FORCE_COMPARISON_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
