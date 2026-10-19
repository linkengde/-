#!/usr/bin/env python3
"""Check full-wavefunction, single-precision GPAW checkpoint write/read."""
import hashlib
import json
import os
import sys
from pathlib import Path

import gpaw
from gpaw import GPAW
from gpaw.mpi import world

mode, source_s, target_s = sys.argv[1:]
source = Path(source_s)
target = Path(target_s)


def sha256(path):
    h = hashlib.sha256()
    fd = os.open(path, os.O_RDONLY)
    try:
        offset = 0
        while True:
            block = os.read(fd, 1024 * 1024)
            if not block:
                break
            h.update(block)
            os.posix_fadvise(fd, offset, len(block), os.POSIX_FADV_DONTNEED)
            offset += len(block)
    finally:
        os.close(fd)
    return h.hexdigest()


if gpaw.__version__ != '26.7.0' or world.size != 4:
    raise RuntimeError(f'unexpected software/ranks: {gpaw.__version__}/{world.size}')

if mode == 'write':
    expected = '38d99727d5b5b3f5152e918685e600790b4a36d819400b552b02519ab92eb2cb'
    if sha256(source) != expected:
        raise RuntimeError('source checkpoint SHA mismatch')
    calc = GPAW(str(source), txt=None)
    calc.dft.ibzwfs.make_sure_wfs_are_read_from_gpw_file()
    if not calc.dft.ibzwfs.has_wave_functions():
        raise RuntimeError('source checkpoint wavefunctions unavailable')
    fd = os.open(source, os.O_RDONLY)
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(fd)
    calc.write(str(target), mode='all', precision='single')
    world.barrier()
    if world.rank == 0:
        print(json.dumps({'stage': 'write', 'status': 'WRITTEN',
                          'source_sha256': expected,
                          'target_bytes': target.stat().st_size,
                          'target_sha256': sha256(target),
                          'precision': 'single', 'mode': 'all',
                          'source_wavefunctions_loaded': True}), flush=True)
elif mode == 'read':
    calc = GPAW(str(target), txt=None)
    calc.dft.ibzwfs.make_sure_wfs_are_read_from_gpw_file()
    if not calc.dft.ibzwfs.has_wave_functions():
        raise RuntimeError('single-precision checkpoint did not restore wavefunctions')
    world.barrier()
    if world.rank == 0:
        print(json.dumps({'stage': 'read', 'status': 'PASS',
                          'target_bytes': target.stat().st_size,
                          'target_sha256': sha256(target),
                          'wavefunctions_reloaded': True}), flush=True)
else:
    raise ValueError(mode)
