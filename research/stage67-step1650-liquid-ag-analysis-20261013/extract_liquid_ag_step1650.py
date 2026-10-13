#!/usr/bin/env python3
"""Reconstruct Stage67 step1650 Liquid-like Ag using archived criteria.

The source dump is read-only. The script writes only a new results directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import ovito
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier, SelectTypeModifier
from scipy.spatial import cKDTree

SOURCE_SHA = "91d04a6fa3cb1b2fd9eb668b392cad04ecb593620723beac3e12ea4e7b8e85fc"
INPUT_COMMIT = "9fa29e3e9dc8c120303f12bcaf8978b701b03604"
REFERENCE_COMMIT = "5fe8b480aaf144ba39e6a9388c8f247e6fadec5a"
EXPECTED_N = 8728
EXPECTED_MAX_CLUSTER = 1032
CLUSTER_CUTOFF = 3.74625
BOX = np.array([150.1, 200.0, 200.0], dtype=float)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def periodic_interval(values: np.ndarray, length: float) -> dict:
    v = np.sort(np.mod(np.asarray(values, dtype=float), length))
    gaps = np.diff(np.r_[v, v[0] + length])
    i = int(np.argmax(gaps))
    start, end = float(v[(i + 1) % len(v)]), float(v[i])
    span = float(length - gaps[i])
    return {
        "start_A": start,
        "end_A": end,
        "span_A": span,
        "center_A": float((start + span / 2) % length),
        "wraps_boundary": bool(start > end and span > 0),
    }


def cluster_candidates(pos: np.ndarray, ids: np.ndarray) -> tuple[list[dict], int]:
    pairs = cKDTree(np.mod(pos, BOX), boxsize=BOX).query_pairs(CLUSTER_CUTOFF, output_type="ndarray")
    parent = np.arange(len(ids), dtype=np.int64)
    weight = np.ones(len(ids), dtype=np.int64)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = int(parent[x])
        return x

    for ai, bi in pairs:
        a, b = find(int(ai)), find(int(bi))
        if a != b:
            if weight[a] < weight[b]:
                a, b = b, a
            parent[b] = a
            weight[a] += weight[b]
    roots = np.fromiter((find(i) for i in range(len(ids))), dtype=np.int64, count=len(ids))
    unique, counts = np.unique(roots, return_counts=True)
    order = np.argsort(counts)[::-1]
    result = []
    for rank, k in enumerate(order):
        m = roots == unique[k]
        result.append({"size": int(counts[k]), "ids": set(map(int, ids[m])), "positions": pos[m], "rank": rank + 1})
    return result, int(len(pairs))


def preserve_filtered_dump(source: Path, target: Path, keep: set[int]) -> int:
    lines = source.read_bytes().splitlines(keepends=True)
    count_i = next(i for i, x in enumerate(lines) if x.strip() == b"ITEM: NUMBER OF ATOMS") + 1
    atoms_i = next(i for i, x in enumerate(lines) if x.startswith(b"ITEM: ATOMS "))
    columns = lines[atoms_i].decode("ascii").split()[2:]
    id_i = columns.index("id")
    records, seen = [], set()
    for line in lines[atoms_i + 1:]:
        fields = line.split()
        if not fields:
            continue
        atom_id = int(fields[id_i])
        if atom_id in seen:
            raise ValueError(f"duplicate source atom ID {atom_id}")
        seen.add(atom_id)
        if atom_id in keep:
            records.append(line)
    if len(seen) != 303372:
        raise ValueError(f"source atom record count mismatch: {len(seen)}")
    ending = b"\r\n" if lines[count_i].endswith(b"\r\n") else b"\n" if lines[count_i].endswith(b"\n") else b""
    lines[count_i] = str(len(records)).encode("ascii") + ending
    target.write_bytes(b"".join(lines[:atoms_i + 1] + records))
    return len(records)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--main-cluster-ids", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    source, summary_path, known_ids_path = (p.resolve(strict=True) for p in (args.source, args.summary, args.main_cluster_ids))
    if args.output_dir.exists():
        raise SystemExit(f"refusing to overwrite: {args.output_dir}")
    if sha256(source) != SOURCE_SHA:
        raise SystemExit("source dump SHA256 mismatch")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    q6 = summary["rows"][0]
    if (q6.get("step"), q6.get("time_ps"), int(q6.get("N_LiquidLikeAg", -1))) != (1650, 0.165, EXPECTED_N):
        raise SystemExit("archived summary does not bind the expected frame/count")

    started = time.time()
    pipe = import_file(str(source), multiple_frames=False)
    pipe.modifiers.append(SelectTypeModifier(types={4}))
    ptm = PolyhedralTemplateMatchingModifier()
    ptm.only_selected = True
    pipe.modifiers.append(ptm)
    data = pipe.compute()
    p = data.particles
    ids = np.asarray(p["Particle Identifier"], dtype=np.int64)
    types = np.asarray(p["Particle Type"], dtype=np.int64)
    selected = np.asarray(p["Selection"], dtype=bool)
    stype = np.asarray(p["Structure Type"], dtype=np.int64)
    pos = np.asarray(p["Position"], dtype=float)
    celltl = np.asarray(p["celltl"], dtype=float)
    valid = np.asarray(p["celltl_valid"], dtype=np.int64)
    ag = types == 4
    thermal_mask = ag & (valid == 1) & (celltl >= 975.0)
    mask = thermal_mask & selected & (stype == 0)
    candidate_ids, candidate_pos = ids[mask], pos[mask]
    if len(candidate_ids) != EXPECTED_N or len(set(map(int, candidate_ids))) != EXPECTED_N:
        raise SystemExit(f"reconstructed candidate count/uniqueness failed: {len(candidate_ids)}")

    clusters, edge_count = cluster_candidates(candidate_pos, candidate_ids)
    known = {int(x) for x in known_ids_path.read_text().splitlines() if x.strip() and not x.lstrip().startswith("#")}
    largest = clusters[0]
    if largest["size"] != EXPECTED_MAX_CLUSTER or largest["ids"] != known:
        raise SystemExit("largest cluster does not exactly match archived 1032-ID list")

    out = args.output_dir.resolve()
    out.mkdir(parents=True)
    out_dump = out / "liquid_Ag_step1650.dump"
    n_written = preserve_filtered_dump(source, out_dump, set(map(int, candidate_ids)))
    if n_written != EXPECTED_N:
        raise SystemExit(f"filtered dump has {n_written} rows, expected {EXPECTED_N}")
    out_ids = out / "liquid_Ag_ids.txt"
    out_ids.write_text("".join(f"{i}\n" for i in sorted(map(int, candidate_ids))), encoding="ascii")

    maxpos = largest["positions"]
    mins, maxs = np.min(maxpos, axis=0), np.max(maxpos, axis=0)
    extents = maxs - mins
    face_gaps = np.minimum(mins, BOX - maxs)
    candidates_envelope = {axis: periodic_interval(candidate_pos[:, i], float(BOX[i])) for i, axis in enumerate("xyz")}
    ptm_other_hot = int((thermal_mask & selected & (stype == 0)).sum())
    ptm_fcc_hot = int((thermal_mask & selected & (stype == 1)).sum())
    ptm_other_ag = int((ag & selected & (stype == 0)).sum())
    ptm_fcc_ag = int((ag & selected & (stype == 1)).sum())
    cluster_sizes = [c["size"] for c in clusters]
    size_json = {
        "status": "PASS",
        "task": "Stage67 step1650 offline Liquid-like Ag extraction",
        "source": {
            "input_commit": INPUT_COMMIT,
            "prior_cluster_commit": REFERENCE_COMMIT,
            "dump_sha256": SOURCE_SHA,
            "dump_bytes": source.stat().st_size,
            "timestep": 1650,
            "time_ps": 0.165,
            "atom_count": 303372,
            "box_bounds_A": {"x": [0.0, 150.1], "y": [0.0, 200.0], "z": [0.0, 200.0]},
            "pbc": [True, True, True],
        },
        "criteria": {
            "description": "Reconstructed from available CellTl/PTM-related criteria in the archived input package.",
            "mask": "Particle Type == 4 AND PTM Structure Type == 0 (Other) AND celltl_valid == 1 AND celltl >= 975 K",
            "ptm_implementation": "OVITO PolyhedralTemplateMatchingModifier",
            "ptm_version": ovito.version_string,
            "ptm_only_selected": True,
            "ptm_rmsd_cutoff": float(ptm.rmsd_cutoff),
            "ptm_enabled_templates": [s.name for s in ptm.structures if s.enabled],
            "q6_threshold_in_summary": q6.get("Q6_threshold"),
            "q6_threshold_applied_to_mask": False,
            "cluster_cutoff_A": CLUSTER_CUTOFF,
        },
        "counts": {
            "all_Ag_atoms": int(ag.sum()),
            "Ag_with_valid_CellTl_ge_975K": int(thermal_mask.sum()),
            "Ag_PTMOther": ptm_other_ag,
            "Ag_PTMFCC": ptm_fcc_ag,
            "hot_valid_Ag_PTMOther_candidates": ptm_other_hot,
            "hot_valid_Ag_PTMFCC_excluded": ptm_fcc_hot,
            "liquid_like_Ag_reconstructed": int(len(candidate_ids)),
            "summary_count_match": int(len(candidate_ids)) == int(q6["N_LiquidLikeAg"]),
            "connected_components": len(clusters),
            "neighbor_edges": edge_count,
        },
        "candidate_region_periodic_envelope": {
            "meaning": "coordinate envelope of all reconstructed candidates, including disconnected islands; not the dimensions of one connected cluster",
            "axis_intervals": candidates_envelope,
        },
        "maximum_cluster": {
            "atom_count": largest["size"],
            "id_list_exact_match_to_archived_maximum_cluster": largest["ids"] == known,
            "center_of_mass_A": np.mean(maxpos, axis=0).tolist(),
            "bounds_min_A": mins.tolist(),
            "bounds_max_A": maxs.tolist(),
            "dimensions_xyz_A": extents.tolist(),
            "distance_to_nearest_simulation_box_face_xyz_A": face_gaps.tolist(),
            "box_faces_are_physical_surfaces": False,
            "PBC_contact": "NO (archived audited value; ID list independently matched)",
            "PBC_spanning": "NO (archived audited value; ID list independently matched)",
        },
        "cluster_statistics": {
            "cluster_count": len(clusters),
            "largest_20_sizes": cluster_sizes[:20],
            "single_atom_cluster_count": sum(n == 1 for n in cluster_sizes),
            "all_sizes_descending": cluster_sizes,
        },
        "validation": {
            "source_hash_pass": True,
            "source_atom_count_pass": True,
            "step_and_time_pass": True,
            "CellTl_fields_present": True,
            "candidate_count_matches_8728_summary": True,
            "maximum_cluster_1032_IDs_exact_match": True,
            "source_atom_records_preserved_byte_for_byte": True,
            "original_dump_modified": False,
            "MD_or_LAMMPS_run": False,
            "offline_analysis_seconds": round(time.time() - started, 3),
        },
    }
    size_path = out / "meltpool_size.json"
    size_path.write_text(json.dumps(size_json, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")

    report = f"""# Stage67 step1650 stable meltpool: final liquid Ag extraction report

## Result

- `FULL_LIQUID_LIKE_AG_EXTRACTION = PASS`
- Reconstructed Liquid-like Ag candidates: **{len(candidate_ids):,}**
- Largest connected candidate cluster: **{largest['size']:,} atoms**
- Connected components at {CLUSTER_CUTOFF:.5f} Å: **{len(clusters)}**
- MD/LAMMPS execution: **none**
- Original dump modified: **no**

**Liquid-like Ag classification was reconstructed from available CellTl/PTM-related criteria in the archived input package.**

The 8,728 count is a reconstructed candidate population, not an independent per-atom label file supplied by the archive. The archive lacks original per-atom PTM/Q6 labels and a complete LiquidCandidate ID table. The classification here was regenerated from the original step1650 dump with the archived reproducible mask; it was not inferred from the 1,032-atom cluster.

## Source and method

- Input package commit: `{INPUT_COMMIT}`; prior maximum-cluster report commit: `{REFERENCE_COMMIT}`.
- Frame: step 1650, t = 0.165 ps; source dump SHA-256 `{SOURCE_SHA}`.
- Model: Stage67 Ag–Ti₃SiC₂ transverse model; 303,372 atoms; cell 150.1 × 200 × 200 Å; `pp pp pp` periodicity.
- Ag type: 4. The dump includes `CellTl` and `CellTl_valid` per atom.
- Reconstructed mask: Ag type 4, OVITO PTM `Structure Type == 0` (`Other`), `CellTl_valid == 1`, and `CellTl >= 975 K`.
- PTM implementation: OVITO {ovito.version_string}, RMSD cutoff {ptm.rmsd_cutoff}, selected Ag centers only. The archived Q6 summary records Q6 threshold {q6.get('Q6_threshold')}; its accompanying extraction report states this value was summary metadata and the per-atom candidate mask used PTM Other. No threshold was changed or added.

The CellTl-valid/hot subset contains {int(thermal_mask.sum()):,} Ag atoms: {ptm_fcc_hot} are PTM FCC and excluded, and {ptm_other_hot:,} are PTM Other candidates. This exactly reproduces the archived summary count of 8,728.

## Cluster and dimensions

PBC-aware connectivity at the archived cutoff {CLUSTER_CUTOFF:.5f} Å produces {len(clusters)} components. The maximum component has {largest['size']:,} Ag atoms, and its complete ID set exactly matches the archived 1,032-ID list. Its center of mass is ({np.mean(maxpos, axis=0)[0]:.6f}, {np.mean(maxpos, axis=0)[1]:.6f}, {np.mean(maxpos, axis=0)[2]:.6f}) Å; its X/Y/Z dimensions are ({extents[0]:.6f}, {extents[1]:.6f}, {extents[2]:.6f}) Å. Distances to the nearest computational-box faces are ({face_gaps[0]:.6f}, {face_gaps[1]:.6f}, {face_gaps[2]:.6f}) Å; these are not physical free surfaces because all three directions are periodic. The audited maximum-cluster flags are `PBC_contact=NO` and `PBC_spanning=NO`.

The `meltpool_size.json` also records a minimum periodic coordinate envelope for all 8,728 candidates. This aggregate contains disconnected components and must not be read as the size of one connected melt pool.

## Integrity and scope

All four uploaded ZIP parts, the reconstructed ZIP, and the package SHA256 list passed. The original dump SHA, 303,372 atom count, step, and CellTl fields passed validation. The output dump preserves selected source atom records byte-for-byte and in source order, changing only the atom-count header to 8,728. IDs are listed in ascending order. The data represent the archived step1650 frame and supplied 0.155–0.165 ps stability window; the package does not include a full time-resolved trajectory or original per-atom classification table. No longer-time or thermodynamic stability claim is made.

See `MANIFEST.json` and `SHA256SUMS.txt` for artifact hashes and analysis environment.
"""
    report_path = out / "STABLE_MELTPOOL_FINAL_REPORT.md"
    report_path.write_text(report, encoding="utf-8")

    script = Path(__file__).resolve()
    manifest = {
        "task": "Stage67 step1650 reconstructed Liquid-like Ag extraction",
        "input_commit": INPUT_COMMIT,
        "reference_commit": REFERENCE_COMMIT,
        "source_dump_sha256": SOURCE_SHA,
        "source_dump_bytes": source.stat().st_size,
        "python": sys.version.split()[0],
        "ovito": ovito.version_string,
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "extractor_sha256": sha256(script),
        "reconstructed_candidate_count": len(candidate_ids),
        "max_cluster_count": largest["size"],
        "max_cluster_ids_match": True,
        "outputs": {
            p.name: sha256(p)
            for p in (out_dump, out_ids, size_path, report_path)
        },
        "original_dump_modified": False,
        "MD_run": False,
    }
    (out / "MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    checksum_files = sorted(p for p in out.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt")
    (out / "SHA256SUMS.txt").write_text("".join(f"{sha256(p)}  {p.name}\n" for p in checksum_files), encoding="ascii")
    print(json.dumps({"status": "PASS", "output_dir": str(out), "candidate_count": len(candidate_ids), "cluster_count": len(clusters), "largest_cluster": largest["size"], "ids_match": True}, indent=2))


if __name__ == "__main__":
    main()
