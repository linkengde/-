#!/usr/bin/env python3
"""Locked Ti3SiC2 48-atom GPAW segmented SCF, dry-run unless --execute.

start/resume/verify/stop only.  No DFT is run by verify or stop.  A source
SHA and approved-H-config SHA are pinned.  The checkpoint scheme requires
actual GPAW validation on H before production; RESOURCE_READY defaults NO.
"""
from __future__ import annotations
import argparse
from collections import Counter
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

INPUT_SHA = 'da3820c477023e0bb423f87550ed82ae0c2a571f53d294c055f78d7d052d5e8b'
CONFIG_SHA = 'ca6f365e8b0b8ec847e329d0312ebb51381b21e996050949494f2d3e0669dfe7'
PAW_SHA = {'Ti': 'bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f',
           'Si': '2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0',
           'C': 'fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983'}
SCOPE = 'TSC48_STATIC_SINGLEPOINT'
MESH = (10, 10, 4)

class Blocked(Exception):
    pass

class Pause(Exception):
    pass

def need(ok, message):
    if not ok:
        raise Blocked(message)

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()

def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def write_atomic(path, data):
    path = Path(path)
    temp = path.with_name(path.name + '.pending')
    need(not temp.exists(), 'STALE_TEMP_FILE ' + str(temp))
    with temp.open('x', encoding='utf-8') as f:
        json.dump(data, f, indent=2, sort_keys=True, ensure_ascii=False)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    os.replace(temp, path)
    sync_directory(path.parent)

def sync_directory(folder):
    if os.name == 'posix':
        fd = os.open(str(folder), os.O_RDONLY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)

def validate_input(source, config):
    need(source.is_file() and config.is_file(), 'INPUT_OR_H_CONFIG_MISSING')
    need(sha256(source) == INPUT_SHA, 'INPUT_SHA_MISMATCH')
    need(sha256(config) == CONFIG_SHA, 'H_CONFIG_SHA_MISMATCH')
    lines = source.read_text(encoding='utf-8').splitlines()
    need(len(lines) == 50 and lines[0].strip() == '48', 'NOT_48_ATOM_EXTXYZ')
    need('pbc="T T T"' in lines[1], 'PBC_MISMATCH')
    need(Counter(line.split()[0] for line in lines[2:]) ==
         Counter({'Ti':24, 'Si':8, 'C':16}), 'COMPOSITION_MISMATCH')
    p = read_json(config)['calculator']
    k = p['kpoints']
    need(k['size'] == [10,10,4] and k['gamma'] == 'omitted / default None'
         and k['offset'] == [0,0,0], 'MP_CONVENTION_MISMATCH')
    need(p['xc'] == 'PBE' and p['mode'] == {'name':'PW','cutoff_eV':500}
         and p['occupations']['width_eV'] == .10
         and p['spin_polarized'] is False, 'METHOD_MISMATCH')
    need(p['eigensolver']['parameters']['niter'] == 2
         and p['mixer']['parameters']['beta'] == .05
         and p['mixer']['parameters']['nmaxold'] == 8
         and p['mixer']['parameters']['weight'] == 100
         and p['maxiter'] == 160
         and p['convergence'] == {'energy':1e-5,'density':1e-5,'eigenstates':1e-5},
         'SCF_MISMATCH')
    need(p['paw_file_sha256'] == {s+'.PBE.gz': h for s,h in PAW_SHA.items()},
         'PAW_CONFIG_MISMATCH')
    return {'scope': SCOPE, 'input_sha256': INPUT_SHA,
            'h_config_sha256': CONFIG_SHA, 'runner_sha256': sha256(__file__)}

def verified_slots(directory, lock):
    valid = []
    for slot in ('A', 'B'):
        meta_file, gpw = directory / ('slot_'+slot+'.json'), directory / ('slot_'+slot+'.gpw')
        if not (meta_file.is_file() and gpw.is_file()):
            continue
        try:
            meta = read_json(meta_file)
            if any(meta.get(k) != v for k,v in lock.items()):
                continue
            if meta['slot'] != slot or meta['mode'] != 'all':
                continue
            if gpw.stat().st_size != meta['bytes'] or sha256(gpw) != meta['sha256']:
                continue
            valid.append((meta, gpw))
        except (ValueError, OSError, KeyError, TypeError):
            pass
    return sorted(valid, key=lambda item: item[0]['sequence'], reverse=True)

def approved_resources(path, directory, seconds, reserve):
    need(path is not None and path.is_file(), 'RESOURCE_CLEARANCE_MISSING')
    c = read_json(path)
    need(c.get('scope') == SCOPE and c.get('resource_ready') is True, 'RESOURCE_READY_NO')
    for name in ('persistent_storage_verified','checkpoint_recovery_validated',
                 'scratch_budget_verified','wallclock_verified'):
        need(c.get(name) is True, 'NOT_VERIFIED_'+name)
    for name in ('persistence_evidence','recovery_evidence','scratch_evidence','wallclock_evidence'):
        need(isinstance(c.get(name),str) and bool(c[name].strip()), 'MISSING_'+name)
    need(Path(c['work_root']).resolve() == directory.parent, 'UNAPPROVED_WORK_ROOT')
    for name in ('wallclock_limit_seconds','checkpoint_write_reserve_seconds','minimum_free_bytes'):
        need(type(c.get(name)) is int and c[name] > 0, 'MISSING_RESOURCE_NUMBER_'+name)
    need(seconds+reserve+c['checkpoint_write_reserve_seconds'] <= c['wallclock_limit_seconds'],
         'NO_VERIFIED_WALLCLOCK_MARGIN')
    need(shutil.disk_usage(directory.parent).free >= c['minimum_free_bytes'],
         'INSUFFICIENT_SITE_CLEARED_FREE_SPACE')
    return c

def root_action(world, directory, callback):
    error = directory / 'root_operation_error.txt'
    if world.rank == 0:
        if error.exists():
            error.unlink()
        try:
            callback()
        except Exception as exc:
            error.write_text(repr(exc), encoding='utf-8')
    world.barrier()
    need(not error.exists(), 'ROOT_OPERATION_FAILED: '+(error.read_text() if error.exists() else ''))

def save_checkpoint(calc, world, directory, lock, sequence, iteration):
    # All ranks enter the GPAW collective write.  No rank-0-only calc.write().
    slot = ('A','B')[sequence % 2]
    stage = directory / ('slot_'+slot+'.writing.gpw')
    final = directory / ('slot_'+slot+'.gpw')
    need(not stage.exists(), 'STALE_STAGE_FILE')
    world.barrier()
    calc.write(str(stage), mode='all')
    world.barrier()
    def publish():
        need(stage.is_file() and stage.stat().st_size > 0, 'CHECKPOINT_WRITE_MISSING')
        with stage.open('rb') as f:
            os.fsync(f.fileno())
        info = {**lock, 'slot':slot, 'sequence':sequence, 'mode':'all',
                'bytes':stage.stat().st_size, 'sha256':sha256(stage),
                'iteration_in_segment':iteration}
        # Only the older slot is overwritten; most recent good slot is retained.
        os.replace(stage, final)
        sync_directory(directory)
        write_atomic(directory / ('slot_'+slot+'.json'), info)
        write_atomic(directory/'segment_status.json',
                     {'status':'CHECKPOINT_WRITTEN', 'sequence':sequence,
                      'checkpoint_sha256':info['sha256']})
    root_action(world, directory, publish)

def run_segment(args, lock):
    import numpy as np
    import gpaw, gpaw_data, ase
    from ase.io import read
    from gpaw import GPAW, PW, restart
    from gpaw.dft import PPCG
    from gpaw.mixer import Mixer
    from gpaw.occupations import FermiDirac
    from gpaw.mpi import world
    need(gpaw.__version__ == '26.7.0' and ase.__version__ == '3.29.0',
         'GPAW_ASE_VERSION_MISMATCH')
    need(world.size == 4, 'MPI_MUST_HAVE_FOUR_RANKS')
    root = Path(gpaw_data.datapath()).resolve()
    need(Path(os.environ.get('GPAW_SETUP_PATH','/none')).resolve() == root,
         'GPAW_SETUP_PATH_MISMATCH')
    for sym, expected in PAW_SHA.items():
        need(sha256(root / (sym+'.PBE.gz')) == expected, 'PAW_HASH_MISMATCH_'+sym)
    directory = args.workdir.resolve()
    clearance = approved_resources(args.clearance, directory, args.segment_seconds,
                                   args.margin_seconds)
    if args.command == 'start':
        need(not directory.exists(), 'WORKDIR_ALREADY_EXISTS')
        if world.rank == 0:
            directory.mkdir(parents=False, exist_ok=False)
        world.barrier()
        root_action(world, directory, lambda: write_atomic(directory/'lock.json', lock))
        sequence, checkpoint = 0, None
    else:
        need(directory.is_dir() and read_json(directory/'lock.json') == lock,
             'RESUME_LOCK_MISMATCH')
        need(not (directory/'RESULT.json').exists(), 'ALREADY_CONVERGED')
        good = verified_slots(directory, lock)
        need(bool(good), 'NO_VERIFIED_FULL_WAVEFUNCTION_CHECKPOINT')
        sequence, checkpoint = good[0][0]['sequence'], good[0][1]
    frozen = read(str(args.input), format='extxyz')
    need(len(frozen)==48 and np.all(frozen.pbc), 'ASE_INPUT_MISMATCH')
    if checkpoint is None:
        atoms = frozen
        calc = GPAW(mode=PW(500), xc='PBE', kpts=MESH,
                    occupations=FermiDirac(.10), eigensolver=PPCG(niter=2),
                    mixer=Mixer(.05,8,100),
                    convergence={'energy':1e-5,'density':1e-5,'eigenstates':1e-5},
                    maxiter=160, spinpol=False, txt=str(directory/'gpaw.start.log'))
        atoms.calc = calc
    else:
        atoms, calc = restart(str(checkpoint),
                              txt=str(directory / ('gpaw.resume.%04d.log' % (sequence+1))))
        need(atoms.get_chemical_symbols() == frozen.get_chemical_symbols()
             and np.array_equal(atoms.positions, frozen.positions)
             and np.array_equal(atoms.cell.array, frozen.cell.array)
             and np.array_equal(atoms.pbc, frozen.pbc), 'RESTART_GEOMETRY_MISMATCH')
        p = calc.params
        need(tuple(p.kpts.todict().get('size', ())) == MESH
             and 'gamma' not in p.kpts.todict() and p.mode.name == 'pw'
             and float(p.mode.ecut) == 500 and p.xc.name == 'PBE'
             and p.maxiter == 160, 'RESTART_SCIENTIFIC_PARAMS_MISMATCH')
        need(p.eigensolver.todict().get('niter') == 2, 'RESTART_SOLVER_MISMATCH')
        atoms.calc = calc
    start = time.monotonic()
    local_iter, ckpt_seq = [0], [sequence]
    stop_file = directory/'STOP_REQUESTED'
    # A stop request remains active until a new segment begins.
    if world.rank == 0 and stop_file.exists():
        stop_file.unlink()
    world.barrier()
    def on_step(ctx):
        local_iter[0] += 1
        due = local_iter[0] % args.checkpoint_every == 0
        end = local_iter[0] >= args.segment_iters
        wall = time.monotonic()-start >= (args.segment_seconds-args.margin_seconds-
                                           clearance['checkpoint_write_reserve_seconds'])
        should_stop = end or wall or stop_file.exists()
        if due or should_stop:
            need(shutil.disk_usage(directory).free >= clearance['minimum_free_bytes'],
                 'CHECKPOINT_DISK_RESERVE_DEPLETED')
            ckpt_seq[0] += 1
            save_checkpoint(calc, world, directory, lock, ckpt_seq[0], int(ctx.niter))
        if should_stop:
            raise Pause()
    calc.hooks['scf_step'] = on_step
    try:
        native = float(atoms.get_potential_energy(force_consistent=False))
    except Pause:
        root_action(world, directory, lambda: write_atomic(directory/'segment_status.json',
                     {'status':'PAUSED_WITH_CHECKPOINT','sequence':ckpt_seq[0],
                      'iterations_this_segment':local_iter[0]}))
        if world.rank == 0:
            print(json.dumps({'status':'PAUSED','checkpoint_sequence':ckpt_seq[0]}))
        return
    free = float(atoms.get_potential_energy(force_consistent=True))
    forces = np.asarray(atoms.get_forces(),dtype=float)
    need(np.isfinite([native,free]).all() and np.isfinite(forces).all(),
         'NONFINITE_FINAL_LABEL')
    try:
        stress = np.asarray(atoms.get_stress(voigt=True),dtype=float).tolist()
        stress_state = 'AVAILABLE'
    except Exception as exc:
        stress,stress_state = None,'UNAVAILABLE_'+type(exc).__name__
    payload = {**lock,'status':'SCF_CONVERGED_KPOINT_UNVALIDATED',
               'data_role':'HOLD','atom_count':48,'kmesh':list(MESH),
               'energy_native_eV':native,'free_energy_eV':free,
               'forces_eV_A':forces.tolist(),'stress_voigt_eV_A3':stress,
               'stress_status':stress_state,
               'stress_convention':'ASE Voigt xx yy zz yz xz xy',
               'iterations_this_segment':local_iter[0]}
    root_action(world, directory, lambda: write_atomic(directory/'RESULT.json',payload))
    if world.rank == 0:
        print(json.dumps({'status':payload['status'],'path':str(directory/'RESULT.json')}))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('verify','start','resume','stop'))
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--h-config',type=Path,required=True)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--clearance',type=Path)
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--segment-iters',type=int,default=20)
    parser.add_argument('--checkpoint-every',type=int,default=20)
    parser.add_argument('--segment-seconds',type=int,default=3600)
    parser.add_argument('--margin-seconds',type=int,default=600)
    args=parser.parse_args()
    try:
        need(args.segment_iters>0 and args.checkpoint_every>0
             and args.segment_seconds>0 and args.margin_seconds>0,'INVALID_SEGMENT_LIMITS')
        lock = validate_input(args.input.resolve(),args.h_config.resolve())
        folder = args.workdir.resolve()
        if args.command == 'verify':
            record={**lock,'status':'OFFLINE_INPUT_CONFIG_PASS'}
            if folder.exists():
                need(folder.is_dir() and read_json(folder/'lock.json')==lock,
                     'WORKDIR_LOCK_MISMATCH')
                record['verified_checkpoint_sequences'] = [m['sequence'] for m,_ in
                                                         verified_slots(folder,lock)]
                record['result_present']=(folder/'RESULT.json').exists()
            print(json.dumps(record,indent=2))
            return
        if args.command == 'stop':
            need(folder.is_dir() and read_json(folder/'lock.json')==lock,
                 'STOP_WORKDIR_MISMATCH')
            write_atomic(folder/'STOP_REQUESTED',{'requested_utc':datetime.datetime.now(
                         datetime.timezone.utc).isoformat()})
            print('STOP_REQUESTED: honored at next SCF callback')
            return
        need(args.execute,'EXECUTION_NOT_REQUESTED_USE_--execute')
        run_segment(args,lock)
    except (Blocked,KeyError,ValueError,OSError) as exc:
        print('BLOCKED: '+str(exc),file=sys.stderr)
        raise SystemExit(2)

if __name__=='__main__':
    main()
