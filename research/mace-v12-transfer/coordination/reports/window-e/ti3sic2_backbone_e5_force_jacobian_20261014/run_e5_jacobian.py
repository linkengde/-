#!/usr/bin/env python3
"""Run frozen-geometry LAMMPS run-0 finite differences for tied TSCBOP blocks.

This is a force-response diagnostic only. It never integrates, minimizes, or
changes coordinates. Every perturbation starts from the same SHA-verified
TSCBOP source and modifies one unordered-pair/radial-field block at a time.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path('/workspace/-')
E4 = ROOT / 'research/mace-v12-transfer/coordination/reports/window-e/ti3sic2_backbone_e4_force_jacobian_20261014'
OUT = ROOT / 'research/mace-v12-transfer/coordination/reports/window-e/ti3sic2_backbone_e5_force_jacobian_20261014'
RUN_ROOT = OUT / 'run0_ucx_self'
POT = ROOT / 'research/mace-v12-transfer/coordination/reports/window-b/current_heating_potential_handoff_20261009/potentials/TSCBOP.tersoff'
LAMMPS = Path('/workspace/.venvs/lammps-run0/bin/lmp')
SOURCE_SHA = '729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7'
FIELDS = {'lambda2': 11, 'B': 12, 'lambda1': 15, 'A': 16}
PARAMS = [
    ('A_Ti-C', 'Ti', 'C', 'A'),
    ('B_Ti-C', 'Ti', 'C', 'B'),
    ('lambda1_Ti-C', 'Ti', 'C', 'lambda1'),
    ('lambda2_Ti-C', 'Ti', 'C', 'lambda2'),
    ('A_Ti-Si', 'Ti', 'Si', 'A'),
    ('B_Ti-Si', 'Ti', 'Si', 'B'),
    ('A_Si-C', 'Si', 'C', 'A'),
    ('B_Si-C', 'Si', 'C', 'B'),
]
TYPE_TO_ELEMENT = {1: 'Ti', 2: 'Si', 3: 'C'}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_manifest(folder: Path) -> int:
    n = 0
    for line in (folder / 'SHA256SUMS.txt').read_text().splitlines():
        digest, rel = line.split(maxsplit=1)
        p = folder / rel.strip().lstrip('*')
        if not p.is_file() or sha(p) != digest:
            raise RuntimeError(f'E4 SHA mismatch: {p}')
        n += 1
    return n


def load_cases():
    label_rows = defaultdict(dict)
    with (E4 / 'dft_force_labels.csv').open(newline='') as f:
        for r in csv.DictReader(f):
            atom_id = int(r['atom_id'])
            label_rows[r['case']][atom_id] = {
                'element': r['element'],
                'xyz': tuple(float(r[k]) for k in ('x_A', 'y_A', 'z_A')),
                'dft': tuple(float(r[k]) for k in ('dft_fx_eV_A', 'dft_fy_eV_A', 'dft_fz_eV_A')),
            }
    e4_manifest = json.loads((E4 / 'dft_label_manifest.json').read_text())
    cases = {}
    for meta in e4_manifest:
        name = meta['case']
        if meta['scf_converged'] is not True or meta['force_unit'] != 'eV/Angstrom':
            raise RuntimeError(f'Invalid DFT force label manifest: {name}')
        atoms = label_rows[name]
        n = int(meta['atom_count'])
        if len(atoms) != n or sorted(atoms) != list(range(1, n + 1)):
            raise RuntimeError(f'Atom ID/count mismatch: {name}')
        if any(not np.isfinite(v) for a in atoms.values() for v in (*a['xyz'], *a['dft'])):
            raise RuntimeError(f'Non-finite label/coordinate: {name}')
        base = E4 / 'run0' / name / 'base'
        data = base / 'atoms.data'
        wrap = base / 'coordinate_wrap_map.csv'
        base_forces = base / 'forces.csv'
        if not all(p.is_file() for p in (data, wrap, base_forces)):
            raise RuntimeError(f'Missing E4 frozen inputs/baseline: {name}')
        # Verify the E4 data file and coordinate map preserve DFT ID/type/position.
        data_atoms = {}
        in_atoms = False
        for line in data.read_text().splitlines():
            if line.strip() == 'Atoms # atomic':
                in_atoms = True
                continue
            if in_atoms and line.strip():
                flds = line.split()
                if len(flds) < 5:
                    continue
                atom_id, typ = int(flds[0]), int(flds[1])
                data_atoms[atom_id] = (TYPE_TO_ELEMENT[typ], tuple(float(x) for x in flds[2:5]))
        if set(data_atoms) != set(atoms):
            raise RuntimeError(f'LAMMPS data ID mapping mismatch: {name}')
        wrapped_by_id = {}
        with wrap.open(newline='') as f:
            for r in csv.DictReader(f):
                i = int(r['atom_id'])
                if i not in atoms or r['element'] != atoms[i]['element']:
                    raise RuntimeError(f'Coordinate map ID/element mismatch: {name}/{i}')
                if any(abs(float(r[k]) - v) > 2e-9 for k, v in zip(('x_original_A','y_original_A','z_original_A'), atoms[i]['xyz'])):
                    raise RuntimeError(f'Coordinate map differs from DFT geometry: {name}/{i}')
                wrapped_by_id[i] = tuple(float(r[k]) for k in ('x_wrapped_A','y_wrapped_A','z_wrapped_A'))
        for i, a in atoms.items():
            el, xyz = data_atoms[i]
            if el != a['element'] or any(abs(x-y) > 2e-9 for x, y in zip(xyz, wrapped_by_id[i])):
                raise RuntimeError(f'Data coordinates differ from DFT geometry: {name}/{i}')
        old = {}
        with base_forces.open(newline='') as f:
            for r in csv.DictReader(f):
                i = int(r['atom_id'])
                old[i] = {'element': r['element'],
                          'xyz': tuple(float(r[k]) for k in ('x_lmp_A','y_lmp_A','z_lmp_A')),
                          'force': tuple(float(r[k]) for k in ('fx_eV_A','fy_eV_A','fz_eV_A'))}
        if set(old) != set(atoms):
            raise RuntimeError(f'E4 baseline ID set mismatch: {name}')
        cases[name] = {'atoms': atoms, 'meta': meta, 'data': data, 'wrap': wrap, 'e4_base': old,
                       'family': 'bulk_parent' if n == 48 else 'interface_derived_slab'}
    if sum(len(c['atoms']) for c in cases.values()) != 354 or len(cases) != 6:
        raise RuntimeError('Unexpected E4 dataset size/case count')
    return cases


def make_variants(source_bytes: bytes):
    text = source_bytes.decode()
    lines = text.splitlines()
    records = []
    for ix, line in enumerate(lines):
        tok = line.split()
        if tok and not tok[0].startswith('#'):
            if len(tok) != 17:
                raise RuntimeError(f'Unexpected Tersoff row width at {ix+1}')
            records.append((ix, tok))
    if len(records) != 27:
        raise RuntimeError(f'Expected 27 Tersoff triplets, got {len(records)}')
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'potentials').mkdir(exist_ok=True)
    source_copy = OUT / 'potentials' / 'TSCBOP.source.tersoff'
    source_copy.write_bytes(source_bytes)
    if sha(source_copy) != SOURCE_SHA:
        raise RuntimeError('Copied source TSCBOP hash mismatch')
    block_rows = []
    variant_records = []
    variants = {}
    for name, el_a, el_b, field in PARAMS:
        idx_field = FIELDS[field]
        pair = {el_a, el_b}
        rows = [(ix, tok) for ix, tok in records if {tok[0], tok[1]} == pair]
        if len(rows) != 6:
            raise RuntimeError(f'{name}: expected six tied pair rows, found {len(rows)}')
        vals = [float(tok[idx_field]) for _, tok in rows]
        if min(vals) <= 0 or max(vals) - min(vals) > 1e-10 * max(vals):
            raise RuntimeError(f'{name}: target field is not one positive shared scalar: {vals}')
        base_value = vals[0]
        variants[name] = {}
        for sign, factor in (('minus1', 0.99), ('plus1', 1.01)):
            out_lines = lines[:]
            for ix, _ in rows:
                t = out_lines[ix].split()
                t[idx_field] = f'{base_value*factor:.12f}'
                out_lines[ix] = ' '.join(t)
            dest = OUT / 'potentials' / f'TSCBOP_{name}_{sign}.tersoff'
            dest.write_text('\n'.join(out_lines) + '\n')
            # Parse-and-compare ensures only the tied scalar column changed.
            changed = []
            for ix, old in records:
                new = out_lines[ix].split()
                if old != new:
                    changed.append((ix, old, new))
            if len(changed) != 6 or any(any(i != idx_field and a[i] != b[i] for i in range(17)) for _, a, b in changed):
                raise RuntimeError(f'{name}/{sign}: variant changed fields beyond its tied block')
            variants[name][sign] = dest
            variant_records.append({'parameter': name, 'sign': sign, 'factor': factor,
                                    'base_value': base_value, 'perturbed_value': base_value*factor,
                                    'path': str(dest.relative_to(OUT)), 'sha256': sha(dest), 'changed_triplet_rows': 6})
        for ix, tok in rows:
            block_rows.append({'parameter': name, 'physical_pair': f'{el_a}-{el_b}', 'field': field,
                               'line_1based': ix+1, 'triplet': ' '.join(tok[:3]), 'base_value': base_value,
                               'minus1_value': base_value*0.99, 'plus1_value': base_value*1.01})
    with (OUT / 'parameter_blocks.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(block_rows[0]), lineterminator='\n'); w.writeheader(); w.writerows(block_rows)
    with (OUT / 'parameter_variants.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(variant_records[0]), lineterminator='\n'); w.writeheader(); w.writerows(variant_records)
    return source_copy, variants, variant_records


def write_input(run_dir: Path, data_path: Path, potential: Path):
    txt = f'''clear
units metal
atom_style atomic
boundary p p p
read_data {data_path.resolve()}
pair_style tersoff
pair_coeff * * {potential.resolve()} Ti Si C
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
thermo 1
thermo_style custom step atoms pe
dump forces all custom 1 forces.dump id type x y z fx fy fz
dump_modify forces sort id format line "%d %d %.16g %.16g %.16g %.16g %.16g %.16g"
run 0 post no
undump forces
'''
    (run_dir / 'run.in').write_text(txt)


def run_case(case_name: str, case: dict, variant: str, potential: Path):
    run_dir = RUN_ROOT / variant.replace('/', '__') / case_name
    run_dir.mkdir(parents=True, exist_ok=True)
    data_path = run_dir / 'atoms.data'
    shutil.copyfile(case['data'], data_path)
    write_input(run_dir, data_path, potential)
    cached = run_dir / 'run_metadata.json'
    if cached.is_file() and (run_dir/'forces.csv').is_file() and (run_dir/'forces.dump').is_file():
        old = json.loads(cached.read_text())
        if (old.get('run_steps') == 0 and old.get('potential_sha256') == sha(potential)
                and old.get('data_sha256') == sha(data_path) and old.get('input_sha256') == sha(run_dir/'run.in')
                and old.get('ucx_tls') == 'self' and old.get('forces_csv_sha256') == sha(run_dir/'forces.csv')
                and old.get('dump_sha256') == sha(run_dir/'forces.dump')):
            force_cache = load_forces_csv(run_dir/'forces.csv')
            if set(force_cache) == set(case['atoms']):
                return {i:{'element':v['element'],'force':v['force']} for i,v in force_cache.items()}, old
    env = os.environ.copy(); env['OMP_NUM_THREADS'] = '1'; env['OMP_DYNAMIC'] = 'FALSE'
    # The singleton MPICH binary defaults to a local UCX socket denied by this
    # sandbox. Self-only transport suffices for one process and leaves the
    # LAMMPS scientific input unchanged.
    env['UCX_TLS'] = 'self'
    cmd = [str(LAMMPS), '-in', 'run.in', '-log', 'lammps.log', '-screen', 'screen.log', '-nocite']
    proc = subprocess.run(cmd, cwd=run_dir, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (run_dir / 'launcher.log').write_text(proc.stdout)
    dump_path = run_dir / 'forces.dump'
    if proc.returncode != 0 or not dump_path.is_file():
        raise RuntimeError(f'LAMMPS run0 failure: {variant}/{case_name}, exit={proc.returncode}')
    lines = dump_path.read_text().splitlines()
    start = next((i for i, x in enumerate(lines) if x.startswith('ITEM: ATOMS')), None)
    if start is None:
        raise RuntimeError(f'Missing force dump atoms section: {variant}/{case_name}')
    names = lines[start].split()[2:]; col = {x:i for i,x in enumerate(names)}
    required = ['id','type','x','y','z','fx','fy','fz']
    if any(k not in col for k in required):
        raise RuntimeError('Incomplete force dump columns')
    forces = {}
    for line in lines[start+1:]:
        if not line.strip():
            continue
        row = line.split(); i = int(row[col['id']]); el = TYPE_TO_ELEMENT[int(row[col['type']])]
        vals = {k: float(row[col[k]]) for k in ('x','y','z','fx','fy','fz')}
        if not all(math.isfinite(v) for v in vals.values()):
            raise RuntimeError(f'Non-finite LAMMPS value: {variant}/{case_name}/{i}')
        forces[i] = {'element': el, 'xyz_lmp': (vals['x'],vals['y'],vals['z']),
                     'force': (vals['fx'],vals['fy'],vals['fz'])}
    atoms = case['atoms']
    if set(forces) != set(atoms):
        raise RuntimeError(f'ID set mismatch: {variant}/{case_name}')
    # LAMMPS wraps the original coordinates into the triclinic box. Compare with E4's
    # audited wrapped coordinates, not a guessed orthogonal-cell transformation.
    wrap_rows = {}
    with case['wrap'].open(newline='') as f:
        for r in csv.DictReader(f):
            wrap_rows[int(r['atom_id'])] = r
    max_coord_err = 0.0
    with (run_dir / 'forces.csv').open('w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['atom_id','element','x_lmp_A','y_lmp_A','z_lmp_A','fx_eV_A','fy_eV_A','fz_eV_A'])
        for i in sorted(atoms):
            a, g = atoms[i], forces[i]
            if g['element'] != a['element']:
                raise RuntimeError(f'Element mismatch {variant}/{case_name}/{i}')
            wr = wrap_rows[i]
            target = tuple(float(wr[k]) for k in ('x_wrapped_A','y_wrapped_A','z_wrapped_A'))
            err = math.sqrt(sum((x-y)**2 for x,y in zip(g['xyz_lmp'], target)))
            max_coord_err = max(max_coord_err, err)
            if err > 1e-8:
                raise RuntimeError(f'PBC coordinate mismatch {variant}/{case_name}/{i}: {err}')
            w.writerow([i,g['element'],*g['xyz_lmp'],*g['force']])
    if max_coord_err > 1e-8:
        raise RuntimeError('Coordinate audit failed')
    log = (run_dir / 'lammps.log').read_text(errors='replace')
    has_step0 = any((parts := line.split()) and parts[0] == '0' and len(parts) >= 3
                    and parts[1] == str(len(atoms)) for line in log.splitlines())
    if 'Loop time' not in log or not has_step0 or 'run 0 post no' not in (run_dir / 'run.in').read_text():
        raise RuntimeError(f'LAMMPS log lacks expected run-0 evidence: {variant}/{case_name}')
    metadata = {
        'case': case_name, 'variant': variant, 'run_steps': 0, 'integrator': False, 'minimizer': False,
        'lammps_version': '22 Jul 2025 Update 4', 'units': 'metal; force=eV/Angstrom',
        'omp_threads': 1, 'ucx_tls': 'self', 'command': cmd, 'potential_sha256': sha(potential),
        'data_sha256': sha(data_path), 'input_sha256': sha(run_dir/'run.in'),
        'dump_sha256': sha(dump_path), 'forces_csv_sha256': sha(run_dir/'forces.csv'),
        'wrapped_coordinate_max_abs_error_A': max_coord_err, 'atom_count': len(atoms),
    }
    (run_dir / 'run_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    return forces, metadata


def load_forces_csv(path: Path):
    out = {}
    with path.open(newline='') as f:
        for r in csv.DictReader(f):
            i = int(r['atom_id'])
            out[i] = {'element': r['element'],
                      'force': tuple(float(r[k]) for k in ('fx_eV_A','fy_eV_A','fz_eV_A'))}
    return out


def main():
    if (OUT / 'run0_status.json').exists():
        raise RuntimeError(f'E5 results already exist; refusing overwrite: {OUT}')
    if not LAMMPS.is_file() or sha(POT) != SOURCE_SHA:
        raise RuntimeError('LAMMPS unavailable or original TSCBOP SHA mismatch')
    e4_count = verify_manifest(E4)
    if json.loads((E4 / 'E4_FINAL_STATUS.json').read_text())['run0_reproduced'] != 'YES':
        raise RuntimeError('E4 baseline reproduction did not pass')
    cases = load_cases()
    source_copy, variants, variant_rows = make_variants(POT.read_bytes())
    (OUT / 'frozen_inputs').mkdir(exist_ok=True)
    for key, case in cases.items():
        target = OUT / 'frozen_inputs' / f'{key}.atoms.data'
        shutil.copyfile(case['data'], target)
        if sha(target) != sha(case['data']):
            raise RuntimeError(f'Frozen data copy hash mismatch: {key}')

    # Recompute one common baseline batch; all parameter columns share these exact
    # geometries, IDs, and LAMMPS settings. First compare against the archived E4 run-0.
    baseline, run_meta = {}, []
    for key, case in cases.items():
        f, m = run_case(key, case, 'baseline', source_copy)
        baseline[key] = f; run_meta.append(m)
        maxdiff = max(abs(f[i]['force'][q] - case['e4_base'][i]['force'][q])
                      for i in f for q in range(3))
        if maxdiff > 1e-8:
            raise RuntimeError(f'E5 source run0 differs from E4 baseline for {key}: {maxdiff}')
        m['max_abs_force_diff_vs_E4_eV_A'] = maxdiff
    (OUT / 'baseline_reproduction.json').write_text(json.dumps({
        'status':'PASS','E4_sha256_manifest_files_verified':e4_count,
        'E4_original_source_sha256':SOURCE_SHA,
        'cases':{k:{'atom_count':len(c['atoms']),'geometry_sha256':c['meta']['geometry_sha256'],
                    'dft_force_source_sha256':c['meta']['dft_force_source_sha256'],
                    'max_abs_force_diff_vs_E4_eV_A':next(m['max_abs_force_diff_vs_E4_eV_A'] for m in run_meta if m['case']==k)}
                 for k,c in cases.items()},
        'force_unit':'eV/Angstrom','run0_only':True}, indent=2)+'\n')

    force_sets = {}
    for name, _, _, _ in PARAMS:
        force_sets[name] = {'minus1':{}, 'plus1':{}}
        for sign in ('minus1','plus1'):
            pot = variants[name][sign]
            for key, case in cases.items():
                f, m = run_case(key, case, f'{name}/{sign}', pot)
                force_sets[name][sign][key] = f
                run_meta.append(m)
    (OUT / 'run0_manifest.json').write_text(json.dumps(run_meta, indent=2)+'\n')

    # The Jacobian rows are stacked in fixed case/ID/component order. Each column
    # differentiates one tied scalar with respect to its log-scaled value.
    atom_response_rows = []
    matrix_rows = []
    per_case_element = []
    per_case = []
    param_arrays = {name: [] for name,_,_,_ in PARAMS}
    residual_dot = defaultdict(float)
    block_sizes = {}
    for key, case in cases.items():
        for i in sorted(case['atoms']):
            a = case['atoms'][i]; b = baseline[key][i]
            vals = {}
            for name,_,_,_ in PARAMS:
                fm = force_sets[name]['minus1'][key][i]['force']
                fp = force_sets[name]['plus1'][key][i]['force']
                j = tuple((fp[q]-fm[q])/0.02 for q in range(3))
                vals[name] = j
                param_arrays[name].extend(j)
                err = tuple(b['force'][q]-a['dft'][q] for q in range(3))
                residual_dot[(name,key,a['element'])] += sum(err[q]*j[q] for q in range(3))
                atom_response_rows.append([name,key,i,a['element'],*a['xyz'],*a['dft'],*b['force'],*fm,*fp,*j,math.sqrt(sum(v*v for v in j))])
            for q, comp in enumerate(('x','y','z')):
                matrix_rows.append([key,i,a['element'],comp,*[vals[name][q] for name,_,_,_ in PARAMS]])
    for name,_,_,_ in PARAMS:
        arr = np.asarray(param_arrays[name], dtype=float)
        block_sizes[name] = {'response_l2_eV_A':float(np.linalg.norm(arr)),
                             'response_rms_component_eV_A':float(np.sqrt(np.mean(arr*arr))),
                             'response_max_abs_component_eV_A':float(np.max(np.abs(arr)))}
    with (OUT/'force_response_by_atom.csv').open('w',newline='') as f:
        cols=['parameter','case','atom_id','element','x_A','y_A','z_A','dft_fx','dft_fy','dft_fz','base_fx','base_fy','base_fz','minus1_fx','minus1_fy','minus1_fz','plus1_fx','plus1_fy','plus1_fz','dFx_dlnp','dFy_dlnp','dFz_dlnp','response_norm']
        w=csv.writer(f,lineterminator='\n');w.writerow(cols);w.writerows(atom_response_rows)
    param_names=[x[0] for x in PARAMS]
    with (OUT/'PARAMETER_RESPONSE_MATRIX.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['case','atom_id','element','component',*param_names]);w.writerows(matrix_rows)

    for name,_,_,_ in PARAMS:
        family_sums = defaultdict(lambda:{'dot':0.0,'response2':0.0,'n':0,'maxnorm':0.0})
        for key,case in cases.items():
            el_acc=defaultdict(lambda:{'j2':0.0,'n':0,'maxnorm':0.0})
            case_j=[]
            for i,a in case['atoms'].items():
                j=tuple((force_sets[name]['plus1'][key][i]['force'][q]-force_sets[name]['minus1'][key][i]['force'][q])/0.02 for q in range(3))
                err=tuple(baseline[key][i]['force'][q]-a['dft'][q] for q in range(3))
                acc=el_acc[a['element']]; acc['j2']+=sum(v*v for v in j); acc['n']+=1; acc['maxnorm']=max(acc['maxnorm'],math.sqrt(sum(v*v for v in j)))
                family=family_sums[(case['family'],a['element'])]
                family['dot']+=sum(err[q]*j[q] for q in range(3)); family['response2']+=sum(v*v for v in j); family['n']+=1; family['maxnorm']=max(family['maxnorm'],math.sqrt(sum(v*v for v in j)))
                case_j.extend(j)
            for el,acc in el_acc.items():
                per_case_element.append([name,key,case['family'],el,acc['n'],math.sqrt(acc['j2']/acc['n']),acc['maxnorm'],residual_dot[(name,key,el)],'increase_reduces_SSE' if residual_dot[(name,key,el)]<0 else ('decrease_reduces_SSE' if residual_dot[(name,key,el)]>0 else 'flat')])
            per_case.append([name,key,case['family'],len(case['atoms']),float(np.sqrt(np.mean(np.square(case_j)))),float(np.linalg.norm(case_j)),max(math.sqrt(sum(v*v for v in ((force_sets[name]['plus1'][key][i]['force'][q]-force_sets[name]['minus1'][key][i]['force'][q])/0.02 for q in range(3)))) for i in case['atoms'] for q in range(3)),sum(residual_dot[(name,key,el)] for el in ('Ti','Si','C'))])
    with (OUT/'force_response_by_element.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['parameter','case','structure_family','element','atom_count','response_rms_vector_eV_A_per_ln_p','max_response_norm_eV_A_per_ln_p','dot_force_residual_J_half_dSSE','local_SSE_direction']);w.writerows(per_case_element)
    with (OUT/'force_response_by_structure.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['parameter','case','structure_family','atom_count','response_rms_component_eV_A_per_ln_p','response_l2_eV_A_per_ln_p','max_abs_component_eV_A_per_ln_p','sum_residual_dot_J_half_dSSE']);w.writerows(per_case)

    J=np.column_stack([np.asarray(param_arrays[n],dtype=float) for n in param_names])
    norms=np.linalg.norm(J,axis=0)
    active_mask=norms > max(norms)*1e-12
    active_names=[n for n,active in zip(param_names,active_mask) if active]
    if not np.any(active_mask):
        raise RuntimeError('All response columns are numerically zero')
    Jn=J[:,active_mask]/norms[active_mask]
    s_raw=np.linalg.svd(J,compute_uv=False)
    s_norm=np.linalg.svd(Jn,compute_uv=False)
    tol_raw=max(J.shape)*np.finfo(float).eps*s_raw[0]
    tol_norm=max(Jn.shape)*np.finfo(float).eps*s_norm[0]
    cosine=Jn.T@Jn
    centered=J-J.mean(axis=0,keepdims=True)
    pnorm=np.linalg.norm(centered,axis=0)
    pearson=np.full((len(param_names),len(param_names)),np.nan)
    for i in range(len(param_names)):
        for j in range(len(param_names)):
            if pnorm[i] > 0 and pnorm[j] > 0:
                pearson[i,j]=float(np.dot(centered[:,i],centered[:,j])/(pnorm[i]*pnorm[j]))
    corr_rows=[]
    for i,a in enumerate(param_names):
        for j,b in enumerate(param_names):
            cos_val=(float(cosine[active_names.index(a),active_names.index(b)]) if a in active_names and b in active_names else None)
            pearson_val=(float(pearson[i,j]) if math.isfinite(pearson[i,j]) else None)
            corr_rows.append([a,b,cos_val,pearson_val,'active' if a in active_names and b in active_names else 'zero_response_column'])
    with (OUT/'parameter_correlation.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['parameter_i','parameter_j','response_cosine','centered_pearson','status']);w.writerows(corr_rows)

    # Compare same-parent perturbations by component-vector cosine; compare the
    # derived slab by per-element response scale and local force-residual direction.
    bulk_cases=[k for k,c in cases.items() if c['family']=='bulk_parent']
    slab_cases=[k for k,c in cases.items() if c['family']=='interface_derived_slab']
    consistency=[]
    for name in param_names:
        for el in ('Ti','Si','C'):
            bulk_vecs={}
            for key in bulk_cases:
                ids=[i for i,a in cases[key]['atoms'].items() if a['element']==el]
                bulk_vecs[key]=np.asarray([(force_sets[name]['plus1'][key][i]['force'][q]-force_sets[name]['minus1'][key][i]['force'][q])/0.02 for i in ids for q in range(3)])
            bulk_cos=[]
            for i,k1 in enumerate(bulk_cases):
                for k2 in bulk_cases[i+1:]:
                    v1,v2=bulk_vecs[k1],bulk_vecs[k2]
                    den=np.linalg.norm(v1)*np.linalg.norm(v2)
                    bulk_cos.append(float(np.dot(v1,v2)/den) if den else 0.0)
            bulk_norm=float(np.sqrt(np.mean(np.concatenate([v*v for v in bulk_vecs.values()]))))
            bulk_dir=np.sign(sum(residual_dot[(name,k,el)] for k in bulk_cases))
            for sk in slab_cases:
                slab_acc=next(r for r in per_case_element if r[0]==name and r[1]==sk and r[3]==el)
                slab_norm=float(slab_acc[5]); slab_dir=np.sign(residual_dot[(name,sk,el)])
                consistency.append([name,el,bulk_norm,slab_norm,slab_norm/(bulk_norm or 1e-300),min(bulk_cos) if bulk_cos else float('nan'),max(bulk_cos) if bulk_cos else float('nan'),bulk_dir,slab_dir,'same' if bulk_dir==slab_dir else ('conflict' if bulk_dir and slab_dir else 'flat')])
    with (OUT/'structure_family_consistency.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['parameter','element','bulk_parent_response_rms_component','derived_slab_response_rms_vector','slab_to_bulk_response_scale','bulk_pairwise_cosine_min','bulk_pairwise_cosine_max','bulk_family_SSE_direction_sign','slab_SSE_direction_sign','direction_comparison']);w.writerows(consistency)

    cond_raw=(float(s_raw[0]/s_raw[-1]) if s_raw[-1] > tol_raw else None)
    cond_norm=(float(s_norm[0]/s_norm[-1]) if s_norm[-1] > tol_norm else None)
    analysis={
        'matrix_shape_components_by_parameters':list(J.shape), 'row_count':int(J.shape[0]), 'column_order':param_names,
        'parameterization':'each column is dF/dln(p), finite difference [F(1.01p)-F(0.99p)]/0.02',
        'raw_jacobian':{'rank_machine_tolerance':int(np.linalg.matrix_rank(J,tol=tol_raw)),
                        'singular_values':s_raw.tolist(),'condition_number':cond_raw,
                        'relative_singular_values':[float(x/s_raw[0]) for x in s_raw],
                        'effective_rank_rel_1e-3':int(np.sum(s_raw/s_raw[0]>1e-3)),
                        'effective_rank_rel_1e-2':int(np.sum(s_raw/s_raw[0]>1e-2)),
                        'machine_rank_tolerance':float(tol_raw)},
        'unit_column_jacobian_active_columns_only':{'column_order':active_names,'rank_machine_tolerance':int(np.linalg.matrix_rank(Jn,tol=tol_norm)),
                                'singular_values':s_norm.tolist(),'condition_number':cond_norm,
                                'relative_singular_values':[float(x/s_norm[0]) for x in s_norm],
                                'effective_rank_rel_1e-3':int(np.sum(s_norm/s_norm[0]>1e-3)),
                                'effective_rank_rel_1e-2':int(np.sum(s_norm/s_norm[0]>1e-2)),
                                'machine_rank_tolerance':float(tol_norm)},
        'column_response_norms':{n:float(norms[i]) for i,n in enumerate(param_names)},
        'numerically_unexcited_columns': [n for n,active in zip(param_names,active_mask) if not active],
        'raw_condition_number_finite':cond_raw, 'active_unit_column_condition_number_finite':cond_norm,
        'data_scope':{'cases':len(cases),'atoms':sum(len(c['atoms']) for c in cases.values()),'bulk_parent_cases':len(bulk_cases),'derived_slab_cases':len(slab_cases),'independent_DFT_family_validation':'not available'},
        'physical_fit_or_optimization_started':False,'DFT_started':False,'MD_or_TTM_started':False,
        'source_tscbop_sha256':SOURCE_SHA,'lammps_version':'22 Jul 2025 Update 4','run0_calls':len(run_meta),
    }
    (OUT/'jacobian_analysis.json').write_text(json.dumps(analysis,indent=2)+'\n')
    (OUT/'run0_status.json').write_text(json.dumps({'status':'PASS','run0_calls':len(run_meta),'baseline_cases':len(cases),'parameter_blocks':len(PARAMS),'variant_runs':len(PARAMS)*2*len(cases),'integrator':False,'minimizer':False,'source_TSCBOP_modified':False,'DFT_started':False,'fit_started':False,'MPI_mode':'single-process LAMMPS binary','OMP_NUM_THREADS':1,'UCX_TLS':'self','note':'self-only transport avoids sandbox-denied local sockets; scientific LAMMPS inputs otherwise identical'},indent=2)+'\n')
    print(json.dumps({'status':'PASS','runs':len(run_meta),'shape':list(J.shape),'rank':analysis['raw_jacobian']['rank_machine_tolerance'],'condition':cond_raw if cond_raw is not None else 'Infinity','active_condition':cond_norm if cond_norm is not None else 'Infinity','active_columns':active_names,'svd_rel':analysis['unit_column_jacobian_active_columns_only']['relative_singular_values'],'columns':analysis['column_response_norms']},indent=2))


if __name__ == '__main__':
    main()
