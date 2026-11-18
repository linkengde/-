#!/usr/bin/env python3
"""Check whether the supplied severe-contact table matches the archived 1M model."""
import csv
import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_IN = Path("/workspace/attachments/531031eb-ec59-44ba-ae1f-1398fdb7185c/full_model_severe_contact_edges.csv")
DATA = Path("/workspace/-/models/ag_ti3sic2_qsgs_1to1/ag_ti3sic2_qsgs.data.gz")
OUT = ROOT / "supplied_severe_contact_csv_assessment.json"
Lx, Ly = 299.635, 301.857484616085
TYPE_ELEMENT = {1: "Ti", 2: "Si", 3: "C", 4: "Ag"}

coords = [None]
types = [0]
in_atoms = False
with gzip.open(DATA, "rt") as f:
    for line in f:
        if line.startswith("Atoms #"):
            in_atoms = True
            continue
        if not in_atoms:
            continue
        p = line.split()
        if len(p) < 5:
            if coords:
                continue
            continue
        try:
            atom_id, typ = int(p[0]), int(p[1])
            xyz = (float(p[2]), float(p[3]), float(p[4]))
        except ValueError:
            if len(coords) > 1:
                break
            continue
        if atom_id >= len(coords):
            coords.extend([None] * (atom_id + 1 - len(coords)))
            types.extend([0] * (atom_id + 1 - len(types)))
        coords[atom_id] = xyz
        types[atom_id] = typ

edges = []
with CSV_IN.open(encoding="utf-8-sig", newline="") as f:
    for row in csv.DictReader(f):
        edges.append((int(row["Ag_ID"]), row["Ag_element"], int(row["TSC_ID"]),
                      row["TSC_element"], float(row["distance_A"])))

in_range = type_match = distance_match = 0
current_distances = []
for aid, ae, bid, be, stated_distance in edges:
    if aid >= len(coords) or bid >= len(coords) or coords[aid] is None or coords[bid] is None:
        continue
    in_range += 1
    if TYPE_ELEMENT.get(types[aid]) != ae or TYPE_ELEMENT.get(types[bid]) != be:
        continue
    type_match += 1
    a, b = coords[aid], coords[bid]
    dx = abs(a[0] - b[0]); dx = min(dx, Lx - dx)
    dy = abs(a[1] - b[1]); dy = min(dy, Ly - dy)
    dz = a[2] - b[2]
    d = (dx * dx + dy * dy + dz * dz) ** 0.5
    current_distances.append((d, stated_distance))
    if abs(d - stated_distance) <= 1e-3:
        distance_match += 1

absdiff = [abs(a - b) for a, b in current_distances]
report = {
    "supplied_csv": str(CSV_IN),
    "comparison_target": str(DATA),
    "supplied_edge_count": len(edges),
    "edge_ids_present_in_target": in_range,
    "edge_ids_with_matching_element_types": type_match,
    "edges_with_distance_match_within_0p001_A": distance_match,
    "median_absolute_distance_difference_A": sorted(absdiff)[len(absdiff)//2] if absdiff else None,
    "supplied_distance_range_A": [min(e[4] for e in edges), max(e[4] for e in edges)],
    "comparison_note": "Target distances use minimum-image x/y periodicity and nonperiodic z, matching the archived static audit boundaries.",
    "classification": "compatible with target" if distance_match >= 0.99 * len(edges) else "not the same atom-id geometry; keep separate until its source structure is identified",
}
OUT.write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
