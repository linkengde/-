#!/usr/bin/env python3
"""Estimate V14 Ag disorder using PTM + Voronoi q6 + displacement."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier, SteinhardtParametersModifier

HERE = Path(__file__).resolve().parent
DEFAULT_START = HERE / "v14_ttm_snapshots/snapshot_000000_start.dump"

def fields(path: Path) -> dict[str, np.ndarray]:
    pipeline = import_file(str(path))
    pipeline.modifiers.append(PolyhedralTemplateMatchingModifier())
    q = SteinhardtParametersModifier()
    q.max_l = 6
    q.neighbor_mode = SteinhardtParametersModifier.NeighborMode.Voronoi
    pipeline.modifiers.append(q)
    p = pipeline.compute().particles
    return {
        "id": np.asarray(p["Particle Identifier"], dtype=np.int64),
        "type": np.asarray(p["Particle Type"], dtype=np.int32),
        "position": np.asarray(p["Position"], dtype=float),
        "ptm": np.asarray(p["Structure Type"], dtype=np.int32),
        "q6": np.asarray(p["Steinhardt"], dtype=float)[:, 4],
    }

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--start", type=Path, default=DEFAULT_START)
    ap.add_argument("--snapshot", type=Path, required=True)
    ap.add_argument("--time-ps", type=float, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--q6-cutoff", type=float, default=0.40)
    ap.add_argument("--msd-cutoff-A", type=float, default=2.0)
    args = ap.parse_args()
    start = fields(args.start)
    if args.start.resolve() == args.snapshot.resolve():
        end = {key: value.copy() for key, value in start.items()}
    else:
        end = fields(args.snapshot)
    for data in (start, end):
        ix = np.argsort(data["id"])
        for k in data:
            data[k] = data[k][ix]
    if not np.array_equal(start["id"], end["id"]):
        raise RuntimeError("Start and snapshot atom ID sets differ")
    p0, p1 = start["position"], end["position"]
    delta = p1 - p0
    delta[:, :2] -= 200.0 * np.rint(delta[:, :2] / 200.0)
    msd2 = np.einsum("ij,ij->i", delta, delta)
    ag = (start["type"] == 2) & (p0[:, 2] >= 40.0) & (p0[:, 2] <= 200.0)
    ptm_other = end["ptm"] == int(PolyhedralTemplateMatchingModifier.Type.OTHER)
    q6_low = end["q6"] < args.q6_cutoff
    mobile = msd2 > args.msd_cutoff_A**2
    candidate = ag & ptm_other & q6_low & mobile
    r = np.hypot((p0[:, 0] - 100.0 + 100.0) % 200.0 - 100.0,
                 (p0[:, 1] - 100.0 + 100.0) % 200.0 - 100.0)
    depth = 200.0 - p0[:, 2]
    profiles = {}
    for label, coord, nbin in (("radial", r, int(200 / np.sqrt(2) / 10) + 1),
                               ("depth", depth, 16)):
        out = []
        bins = np.floor(np.maximum(coord, 0.0) / 10.0).astype(int)
        for b in range(nbin):
            cell = ag & (bins == b)
            n = int(cell.sum())
            if n >= 30:
                bounds = [b*10.0, (b+1)*10.0]
                out.append({"range_A": bounds, "Ag_atoms": n,
                            "melt_candidate_fraction": float((candidate & cell).sum()/n)})
        profiles[label] = out
    ag_n = int(ag.sum())
    report = {
        "status": "EXPLORATORY_STRUCTURAL_ESTIMATE",
        "time_ps": args.time_ps,
        "inputs": {"start": str(args.start), "snapshot": str(args.snapshot)},
        "thresholds": {"PTM": "OTHER", "Voronoi_q6_below": args.q6_cutoff,
                       "displacement_above_A": args.msd_cutoff_A},
        "Ag_initial_TTM_layer": {
            "atoms": ag_n,
            "PTM_OTHER": int((ag & ptm_other).sum()),
            "q6_below_cutoff": int((ag & q6_low).sum()),
            "displacement_above_cutoff": int((ag & mobile).sum()),
            "combined_candidates": int(candidate.sum()),
            "combined_fraction": float(candidate.sum()/ag_n) if ag_n else 0.0,
            "candidate_rms_displacement_A": float(np.sqrt(msd2[candidate].mean())) if candidate.any() else 0.0,
        },
        "profiles_10A": profiles,
        "caveats": [
            "Strong disorder classifier, not an equilibrium liquid coexistence test.",
            "Spatial profiles use initial positions; large atom motion can blur the melt boundary.",
            "TTM map is frozen during the run; phase migration is not remapped dynamically.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["Ag_initial_TTM_layer"], indent=2))

if __name__ == "__main__":
    main()
