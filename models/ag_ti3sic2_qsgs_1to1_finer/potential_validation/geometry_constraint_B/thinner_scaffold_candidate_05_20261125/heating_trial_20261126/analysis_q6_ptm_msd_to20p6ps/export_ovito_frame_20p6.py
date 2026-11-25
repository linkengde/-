#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用OVITO读取20.6 ps帧并导出坐标、元素类型和Ag局部qbar6。"""
from pathlib import Path
import numpy as np
from ovito.io import import_file
from analyze_q6_ptm_msd import read_ag

HERE = Path(__file__).resolve().parent
FRAME = HERE / "frame_20.600ps_step41200.dump.gz"
pipeline = import_file(str(FRAME))
data = pipeline.compute()
ids = np.asarray(data.particles["Particle Identifier"].array, dtype=np.int64)
types = np.asarray(data.particles["Particle Type"].array, dtype=np.int32)
positions = np.asarray(data.particles["Position"].array, dtype=np.float32)
ag_ids, _ag_pos, _q6, qbar6, _ptm, _names = read_ag(FRAME)
qbar_by_id = {int(i): float(q) for i, q in zip(ag_ids, qbar6)}
qbar_all = np.full(len(ids), np.nan, dtype=np.float32)
ag = types == 4
qbar_all[ag] = np.fromiter((qbar_by_id[int(i)] for i in ids[ag]), dtype=np.float32, count=int(ag.sum()))
out = HERE / "frame_20p6ps_ovito_analysis_arrays.npz"
np.savez_compressed(out, positions=positions, types=types, particle_ids=ids, qbar6=qbar_all)
print(f"WROTE {out}: atoms={len(ids)}, Ag={int(ag.sum())}", flush=True)
