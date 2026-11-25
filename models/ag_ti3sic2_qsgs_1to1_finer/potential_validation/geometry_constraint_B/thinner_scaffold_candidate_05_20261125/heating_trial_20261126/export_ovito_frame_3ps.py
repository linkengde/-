#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从 3.0 ps LAMMPS 帧用 OVITO 导出坐标、类型与 Ag-only qbar6。"""
from pathlib import Path
import sys

import numpy as np
from ovito.io import import_file

HERE = Path(__file__).resolve().parent
FRAME = HERE / "continued_2p5_to_3p0ps" / "frame_03.000ps_step06000.dump.gz"
sys.path.insert(0, str(HERE / "analysis_q6_ptm_msd_0_to_3ps"))
from analyze_q6_ptm_msd import read_ag  # noqa: E402

pipeline = import_file(str(FRAME))
data = pipeline.compute()
ids = np.asarray(data.particles["Particle Identifier"].array, dtype=np.int64)
types = np.asarray(data.particles["Particle Type"].array, dtype=np.int32)
positions = np.asarray(data.particles["Position"].array, dtype=np.float32)
ag_ids, _ag_pos, _q6, qbar6, _ptm, _names = read_ag(FRAME)
qbar_by_id = {int(i): float(q) for i, q in zip(ag_ids, qbar6)}
qbar_all = np.full(len(ids), np.nan, dtype=np.float32)
ag = types == 4
qbar_all[ag] = np.fromiter((qbar_by_id[int(i)] for i in ids[ag]),
                           dtype=np.float32, count=int(np.sum(ag)))
out = HERE / "frame_3ps_ovito_analysis_arrays.npz"
np.savez_compressed(out, positions=positions, types=types, particle_ids=ids, qbar6=qbar_all)
print(f"Wrote {out}; N={len(ids)}, Ag={int(np.sum(ag))}, skeleton={int(np.sum(~ag))}", flush=True)
