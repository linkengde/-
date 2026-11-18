#!/workspace/venvs/dft-labels/bin/python
"""Stop the conservative FIRE attempt if any recorded Fmax exceeds 500 eV/A."""
from __future__ import annotations

import json
import os
import signal
import sys
import time
from pathlib import Path

folder = Path(__file__).resolve().parent
pid = int(sys.argv[1])
log = folder / 'log.lammps'
marker = folder / 'watchdog_stop.json'
last_step = -1
while True:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        break
    try:
        lines = log.read_text(errors='replace').splitlines()
    except FileNotFoundError:
        time.sleep(2)
        continue
    for line in lines:
        parts = line.split()
        if len(parts) < 7 or parts[1] != '347464' or not parts[0].isdigit():
            continue
        step = int(parts[0])
        if step <= last_step:
            continue
        last_step = step
        fmax = float(parts[4])
        if fmax > 500.0:
            marker.write_text(json.dumps({'status': 'STOPPED_FORCE_LIMIT', 'pid': pid,
                                          'step': step, 'Fmax_eV_A': fmax,
                                          'signal': 'SIGINT'}, indent=2) + '\n')
            try:
                os.kill(pid, signal.SIGINT)
            except ProcessLookupError:
                pass
            raise SystemExit(0)
    time.sleep(2)
marker.write_text(json.dumps({'status': 'LAMMPS_EXITED_BELOW_FORCE_LIMIT', 'pid': pid,
                              'last_monitored_step': last_step}, indent=2) + '\n')
