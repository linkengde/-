#!/usr/bin/env python3
"""联合计算 Ag 热源区的 q6、邻居平均 qbar6、PTM 和位移 MSD。"""

import csv
import json
import re
from pathlib import Path

import numpy as np
from ovito.io import import_file
from ovito.modifiers import (
    DeleteSelectedModifier,
    PolyhedralTemplateMatchingModifier,
    SelectTypeModifier,
    SteinhardtParametersModifier,
)


HERE = Path(__file__).resolve().parent
GEOMETRY = json.loads((HERE / "geometry_report.json").read_text())
LX, LY, _ = map(float, GEOMETRY["box_A"])
ZTOP = float(GEOMETRY["material_z_A"][1])
CORE_R_A = 30.0
CORE_DEPTH_A = 30.0
Q6_CUTOFF_A = 3.6
Q6_SCREEN_THRESHOLD = 0.45


def frame_time(path):
    """从标准输出帧名中解析时间；step0 初态时间为 0。"""
    if path.name.startswith("frame_00000"):
        return 0.0
    match = re.search(r"frame_(\d+\.\d+)ps_", path.name)
    if not match:
        raise ValueError(f"无法从帧名解析时间：{path.name}")
    return float(match.group(1))


def read_ag(path):
    """删除 Ti/Si/C 后，在纯 Ag 邻居图上算 q6、qbar6 和 PTM。"""
    pipeline = import_file(str(path))
    pipeline.modifiers.append(SelectTypeModifier(types={1, 2, 3}))
    pipeline.modifiers.append(DeleteSelectedModifier())
    qmod = SteinhardtParametersModifier()
    qmod.max_l = 6
    qmod.neighbor_mode = SteinhardtParametersModifier.NeighborMode.Cutoff
    qmod.cutoff = Q6_CUTOFF_A
    qmod.compute_averaged = True
    pipeline.modifiers.append(qmod)
    pipeline.modifiers.append(PolyhedralTemplateMatchingModifier())
    data = pipeline.compute()

    ids = np.asarray(data.particles["Particle Identifier"].array, dtype=np.int64)
    pos = np.asarray(data.particles["Position"].array, dtype=np.float64)
    order = np.argsort(ids)
    q6_prop = data.particles["Steinhardt"]
    q6_index = q6_prop.component_names.index("q6")
    qbar_prop = data.particles["Smooth Steinhardt"]
    qbar_index = qbar_prop.component_names.index("q6")
    q6 = np.asarray(q6_prop.array[:, q6_index], dtype=np.float64)[order]
    qbar6 = np.asarray(qbar_prop.array[:, qbar_index], dtype=np.float64)[order]
    ptm = np.asarray(data.particles["Structure Type"].array, dtype=np.int32)[order]
    names = {int(t.id): t.name for t in data.particles["Structure Type"].types}
    return ids[order], pos[order], q6, qbar6, ptm, names


def order_stats(values):
    """返回连续分布统计和一个仅供相对比较的 q6>0.45 比例。"""
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "p10": float(np.percentile(values, 10)),
        "p90": float(np.percentile(values, 90)),
        "fraction_above_screen_threshold": float(np.mean(values > Q6_SCREEN_THRESHOLD)),
    }


def main():
    files = [p for p in HERE.glob("frame_*.dump.gz") if p.name.startswith("frame_")]
    files.sort(key=frame_time)
    if not files or frame_time(files[0]) != 0.0:
        raise FileNotFoundError("分析目录必须包含 frame_00000_step0.dump.gz")

    records = []
    by_time = {}
    initial_ids = None
    initial_pos = None
    core_ids = None

    for path in files:
        time_ps = frame_time(path)
        ids, pos, q6, qbar6, ptm, names = read_ag(path)
        if initial_ids is None:
            initial_ids = ids.copy()
            initial_pos = pos.copy()
            r2 = (pos[:, 0] - LX / 2.0) ** 2 + (pos[:, 1] - LY / 2.0) ** 2
            core_mask0 = (
                (r2 <= CORE_R_A**2)
                & (pos[:, 2] >= ZTOP - CORE_DEPTH_A)
                & (pos[:, 2] <= ZTOP)
            )
            core_ids = ids[core_mask0]
        if not np.array_equal(ids, initial_ids):
            raise RuntimeError(f"Ag 粒子 ID 改变：{path.name}")

        core_mask = np.isin(ids, core_ids)
        for region_name, mask in (
            ("all_Ag", np.ones(len(ids), dtype=bool)),
            ("XY_core_top30A", core_mask),
        ):
            vals_q6 = q6[mask]
            vals_qbar6 = qbar6[mask]
            vals_ptm = ptm[mask]
            row = {
                "time_ps": time_ps,
                "frame": path.name,
                "region": region_name,
                "n_Ag": int(mask.sum()),
                "q6_mean": float(np.mean(vals_q6)),
                "q6_median": float(np.median(vals_q6)),
                "q6_p10": float(np.percentile(vals_q6, 10)),
                "q6_p90": float(np.percentile(vals_q6, 90)),
                "q6_fraction_gt_0p45": float(np.mean(vals_q6 > Q6_SCREEN_THRESHOLD)),
                "qbar6_mean": float(np.mean(vals_qbar6)),
                "qbar6_median": float(np.median(vals_qbar6)),
                "qbar6_p10": float(np.percentile(vals_qbar6, 10)),
                "qbar6_p90": float(np.percentile(vals_qbar6, 90)),
                "qbar6_fraction_gt_0p45": float(np.mean(vals_qbar6 > Q6_SCREEN_THRESHOLD)),
                "PTM_other_fraction": float(np.mean(vals_ptm == 0)),
                "PTM_FCC_fraction": float(np.mean(vals_ptm == 1)),
                "PTM_HCP_fraction": float(np.mean(vals_ptm == 2)),
            }
            records.append(row)

        # 计算相对于初始帧的 MSD；xy 方向按周期最短像修正，z 为开放边界。
        # 同时报告去除热区整体平移后的 MSD，避免把刚体漂移误当作扩散。
        delta = pos - initial_pos
        delta[:, 0] -= LX * np.rint(delta[:, 0] / LX)
        delta[:, 1] -= LY * np.rint(delta[:, 1] / LY)
        msd = np.einsum("ij,ij->i", delta, delta)
        core_com_drift = np.mean(delta[core_mask], axis=0)
        delta_com_corrected = delta - core_com_drift
        msd_com_corrected = np.einsum(
            "ij,ij->i", delta_com_corrected, delta_com_corrected
        )
        by_time[time_ps] = {
            "frame": path.name,
            "all_Ag_msd_A2": float(np.mean(msd)),
            "hot_core_msd_A2": float(np.mean(msd[core_mask])),
            "hot_core_rms_displacement_A": float(np.sqrt(np.mean(msd[core_mask]))),
            "hot_core_com_drift_A": [float(v) for v in core_com_drift],
            "hot_core_msd_com_corrected_A2": float(
                np.mean(msd_com_corrected[core_mask])
            ),
            "hot_core_rms_com_corrected_A": float(
                np.sqrt(np.mean(msd_com_corrected[core_mask]))
            ),
        }
        print(f"ANALYZED {time_ps:.3f} ps: q6/PTM/MSD, hot_core_Ag={len(core_ids)}", flush=True)

    csv_path = HERE / "q6_ptm_msd_timeseries.csv"
    fields = list(records[0].keys())
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)

    msd_path = HERE / "msd_timeseries.csv"
    with msd_path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "time_ps", "frame", "all_Ag_msd_A2", "hot_core_msd_A2",
                "hot_core_rms_displacement_A", "hot_core_com_drift_A",
                "hot_core_msd_com_corrected_A2", "hot_core_rms_com_corrected_A",
            ],
        )
        writer.writeheader()
        for time_ps in sorted(by_time):
            writer.writerow({"time_ps": time_ps, **by_time[time_ps]})

    summary = {
        "method": "OVITO PTM + Steinhardt q6/qbar6 + initial-origin MSD; Ag-only neighbor graph",
        "q6_neighbor_rule": f"uniform cutoff {Q6_CUTOFF_A:.2f} A",
        "q6_screen_threshold": Q6_SCREEN_THRESHOLD,
        "q6_threshold_note": "0.45 is only a same-model trend screen, not a universal solid/liquid boundary.",
        "hot_region": f"r_xy<={CORE_R_A:g} A from box center and z>={ZTOP-CORE_DEPTH_A:g} A",
        "hot_region_atom_count_initial": int(len(core_ids)),
        "PTM_codes": names,
        "phase_series": records,
        "MSD_series": [{"time_ps": t, **by_time[t]} for t in sorted(by_time)],
        "decision_rule": (
            "Do not infer melting from PTM Other alone. Look for a coherent hot-core q6/qbar6 decrease, "
            "loss of PTM crystalline motifs, and a sustained rising hot-core MSD together."
        ),
    }
    (HERE / "q6_ptm_msd_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n"
    )
    print(f"WROTE {csv_path} and {msd_path}", flush=True)


if __name__ == "__main__":
    main()
