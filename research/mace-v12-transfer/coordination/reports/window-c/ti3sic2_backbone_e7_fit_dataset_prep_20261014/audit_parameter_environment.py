#!/usr/bin/env python3
"""Geometry-only audit of E5 TSCBOP pair-cutoff activation.

Reads the ten imported Stage C extxyz geometries using ASE neighbor lists.
It does not load/evaluate TSCBOP, construct a calculator, run LAMMPS, or fit.
The radial cutoffs below were transcribed from the E5 archived source file
whose SHA-256 is recorded in provenance.json.
"""

from __future__ import annotations

import csv
from pathlib import Path

from ase.io import read
from ase.neighborlist import neighbor_list


ROOT = Path(__file__).resolve().parent
IMPORT = ROOT / "stage_c_import"
GEOMS = IMPORT / "geometry_only_inputs"
PAIR_CUTOFFS = {
    "Ti-C": {"R_A": 2.731200, "D_A": 0.187530,
             "blocks": "A_Ti-C;B_Ti-C;lambda1_Ti-C;lambda2_Ti-C"},
    "Ti-Si": {"R_A": 3.256477, "D_A": 0.297824,
              "blocks": "A_Ti-Si;B_Ti-Si"},
    "Si-C": {"R_A": 2.960379, "D_A": 0.251126,
             "blocks": "A_Si-C;B_Si-C"},
}
ELEMENT_PAIR = {tuple(sorted(k.split("-"))): k for k in PAIR_CUTOFFS}


def main() -> None:
    with (IMPORT / "DFT_STRUCTURE_LIST.csv").open(newline="", encoding="utf-8") as f:
        inventory = list(csv.DictReader(f))
    rows: list[dict] = []
    for item in inventory:
        if item["structure_status"] != "PASS_GEOMETRY_ONLY":
            for pair, p in PAIR_CUTOFFS.items():
                rows.append({
                    "structure_id": item["structure_id"],
                    "category": item["category"],
                    "structure_status": "SOURCE_PENDING_NOT_EXECUTABLE",
                    "dft_label_status": "NO_STRUCTURE_OR_DFT_LABEL",
                    "natoms": "", "pair": pair,
                    "inner_full_cutoff_R_minus_D_A": f"{p['R_A'] - p['D_A']:.6f}",
                    "outer_cutoff_R_plus_D_A": f"{p['R_A'] + p['D_A']:.6f}",
                    "unique_periodic_pairs_inside_outer_cutoff": "NA",
                    "unique_periodic_pairs_in_taper": "NA",
                    "minimum_pair_distance_A": "NA",
                    "maximum_pair_distance_inside_outer_A": "NA",
                    "mapped_parameter_blocks": p["blocks"],
                    "interpretation": "Cannot assess without an independently sourced structure; this is only a reserved role slot.",
                })
            continue

        path = GEOMS / Path(item["structure_file"]).name
        atoms = read(path, format="extxyz")
        cutoff = max(p["R_A"] + p["D_A"] for p in PAIR_CUTOFFS.values()) + 1.0
        ii, jj, dd = neighbor_list("ijd", atoms, cutoff, self_interaction=False)
        pair_distances = {p: [] for p in PAIR_CUTOFFS}
        for i, j, d in zip(ii, jj, dd):
            # ASE returns directed bonds. Keep one direction; periodic image
            # shifts remain distinct, physically distinct pair interactions.
            if int(i) >= int(j):
                continue
            key = ELEMENT_PAIR.get(tuple(sorted((atoms[int(i)].symbol,
                                                 atoms[int(j)].symbol))))
            if key is not None:
                pair_distances[key].append(float(d))

        for pair, p in PAIR_CUTOFFS.items():
            values = pair_distances[pair]
            inner = p["R_A"] - p["D_A"]
            outer = p["R_A"] + p["D_A"]
            active = [d for d in values if d < outer - 1e-10]
            taper = [d for d in active if d > inner + 1e-10]
            rows.append({
                "structure_id": item["structure_id"],
                "category": item["category"],
                "structure_status": "PASS_GEOMETRY_ONLY",
                "dft_label_status": "NO_DFT_ENERGY_FORCE_OR_STRESS_LABEL",
                "natoms": len(atoms), "pair": pair,
                "inner_full_cutoff_R_minus_D_A": f"{inner:.6f}",
                "outer_cutoff_R_plus_D_A": f"{outer:.6f}",
                "unique_periodic_pairs_inside_outer_cutoff": len(active),
                "unique_periodic_pairs_in_taper": len(taper),
                "minimum_pair_distance_A": f"{min(values):.6f}" if values else "NA",
                "maximum_pair_distance_inside_outer_A": f"{max(active):.6f}" if active else "NA",
                "mapped_parameter_blocks": p["blocks"],
                "interpretation": (
                    "Nonzero geometric cutoff coverage only; not a force-response derivative, identifiability test, or DFT fit label."
                    if active else
                    "No pair enters the nonzero Tersoff radial cutoff; mapped blocks receive zero direct radial pair contribution on this geometry."
                ),
            })

    out = ROOT / "PARAMETER_ENVIRONMENT_ACTIVATION.csv"
    fields = list(rows[0].keys())
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} pair/structure rows to {out}")


if __name__ == "__main__":
    main()
