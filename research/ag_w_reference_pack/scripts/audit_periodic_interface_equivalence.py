#!/usr/bin/env python3
"""Audit whether the two periodic interfaces in each EFS cell are equivalent."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from ase.io import write
from scipy.spatial import cKDTree

from relax_and_compare_ag_w_interfaces import (
    CONVERSION_J_M2,
    INTERFACES,
    POTENTIAL,
    build_interface_atoms,
    build_orientation,
    core_region,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE_REPORT = ROOT / "validation" / "relaxed_interface_report.json"
OUTPUT = ROOT / "validation" / "two_interface_equivalence_audit.json"
CASES = {
    "Ag100_W100": {"registry": "top", "gap_A": 2.8},
    "Ag111_W110": {"registry": "shift_mixed", "gap_A": 2.8},
}
DISTANCE_BINS_A = [0.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.3, 4.6]
WINDOW_HALF_WIDTH_A = 5.0


def read_lammps_dump(path: Path) -> dict[str, np.ndarray]:
    lines = path.read_text().splitlines()
    atom_header = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS"))
    fields = lines[atom_header].split()[2:]
    data = np.loadtxt(lines[atom_header + 1 :], dtype=float)
    return {name: data[:, i] for i, name in enumerate(fields)}


def run_minimized_case(lammps: str, name: str, registry: str, gap: float, work: Path) -> dict:
    spec = INTERFACES[name]
    geometry = build_orientation(name)
    atoms, info = build_interface_atoms(geometry, gap, spec["registries"][registry], spec["core_depth_A"])
    run_dir = work / name
    run_dir.mkdir()
    write(
        run_dir / "atoms.data",
        atoms,
        format="lammps-data",
        atom_style="atomic",
        masses=True,
        specorder=["Ag", "W"],
    )
    regions = "".join(core_region(bounds, region) for region, bounds in info["core_regions"].items())
    group_names = list(info["core_regions"])
    core_group = f"group core union {' '.join(group_names)}"
    (run_dir / "in.lmp").write_text(
        f"""units metal
atom_style atomic
boundary p p p
read_data atoms.data
pair_style eam/fs
pair_coeff * * {POTENTIAL} Ag W
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
{regions}{core_group}
fix hold core setforce 0.0 0.0 0.0
compute e_atom all pe/atom
min_style cg
minimize 0.0 1e-5 8000 80000
write_dump all custom relaxed.dump id type x y z c_e_atom modify sort id
variable etotal equal pe
print "RESULT energy=$(v_etotal:%.14g)"
"""
    )
    result = subprocess.run([lammps, "-in", "in.lmp", "-log", "none"], cwd=run_dir, capture_output=True, text=True)
    output = result.stdout + result.stderr
    if result.returncode:
        raise RuntimeError(f"LAMMPS failed for {name}:\n{output[-3000:]}")
    match = re.search(r"RESULT energy=([-+0-9.eE]+)", output)
    if match is None:
        raise RuntimeError(f"LAMMPS omitted the final energy for {name}:\n{output[-2000:]}")
    dump = read_lammps_dump(run_dir / "relaxed.dump")

    ag_max = float(geometry["ag"].positions[:, 2].max())
    ag_min = float(geometry["ag"].positions[:, 2].min())
    w_min = float(geometry["w"].positions[:, 2].min())
    w_max = float(geometry["w"].positions[:, 2].max())
    w_high = ag_max + gap + w_max - w_min
    box_z = float(atoms.cell.lengths()[2])
    centers = {
        "Ag_top_W_bottom": ag_max + gap / 2.0,
        "W_top_Ag_bottom_periodic": w_high + gap / 2.0,
    }
    orientation_report = json.loads(SOURCE_REPORT.read_text())["orientations"][name]
    bulk_epa = {
        "Ag": orientation_report["surface_references"]["Ag"]["bulk_energy_eV_atom"],
        "W": orientation_report["surface_references"]["W"]["bulk_energy_eV_atom"],
    }
    area = orientation_report["area_A2"]
    result = {
        "registry": registry,
        "initial_gap_A": gap,
        "final_total_energy_eV": float(match.group(1)),
        "window_half_width_A": WINDOW_HALF_WIDTH_A,
        "interfaces": {},
    }

    types = dump["type"].astype(int)
    z = dump["z"] % box_z
    for interface_name, center in centers.items():
        dz = np.abs((z - center + box_z / 2.0) % box_z - box_z / 2.0)
        selected = dz < WINDOW_HALF_WIDTH_A
        atom_counts = {
            "Ag": int(np.count_nonzero(selected & (types == 1))),
            "W": int(np.count_nonzero(selected & (types == 2))),
        }
        excess = sum(
            float(np.sum(dump["c_e_atom"][selected & (types == atom_type)] - bulk_epa[element]))
            for atom_type, element in ((1, "Ag"), (2, "W"))
        )
        local = {
            "center_z_A": float(center),
            "atoms_within_window": atom_counts,
            "local_excess_energy_per_area_J_m2": excess / area * CONVERSION_J_M2,
        }
        result["interfaces"][interface_name] = local

    box = np.array([atoms.cell.lengths()[0], atoms.cell.lengths()[1], box_z])
    pos = np.column_stack((dump["x"] % box[0], dump["y"] % box[1], z))
    pos[pos >= box - 1.0e-10] = 0.0
    ag_pos = pos[types == 1]
    w_pos = pos[types == 2]
    pairs = cKDTree(w_pos, boxsize=box).sparse_distance_matrix(
        cKDTree(ag_pos, boxsize=box), 4.6, output_type="ndarray"
    )
    cross_dz = ag_pos[pairs["j"], 2] - w_pos[pairs["i"], 2]
    cross_dz -= np.rint(cross_dz / box_z) * box_z
    mid_z = (w_pos[pairs["i"], 2] + cross_dz / 2.0) % box_z
    for interface_name, center in centers.items():
        center_delta = np.abs((mid_z - center + box_z / 2.0) % box_z - box_z / 2.0)
        in_window = center_delta < WINDOW_HALF_WIDTH_A
        distances = pairs["v"][in_window]
        result["interfaces"][interface_name].update(
            {
                "Ag_W_pairs_within_4p6_A": int(len(distances)),
                "minimum_Ag_W_distance_A": float(distances.min()),
                "Ag_W_pair_distance_histogram": np.histogram(distances, bins=DISTANCE_BINS_A)[0].astype(int).tolist(),
            }
        )

    a = result["interfaces"]["Ag_top_W_bottom"]
    b = result["interfaces"]["W_top_Ag_bottom_periodic"]
    energy_difference = abs(a["local_excess_energy_per_area_J_m2"] - b["local_excess_energy_per_area_J_m2"])
    pair_total = max(a["Ag_W_pairs_within_4p6_A"], b["Ag_W_pairs_within_4p6_A"], 1)
    pair_count_difference = abs(a["Ag_W_pairs_within_4p6_A"] - b["Ag_W_pairs_within_4p6_A"]) / pair_total
    hist_a = np.asarray(a["Ag_W_pair_distance_histogram"], dtype=float)
    hist_b = np.asarray(b["Ag_W_pair_distance_histogram"], dtype=float)
    hist_l1_fraction = float(np.abs(hist_a - hist_b).sum() / max(hist_a.sum() + hist_b.sum(), 1.0))
    result["comparison"] = {
        "absolute_local_excess_energy_difference_J_m2": energy_difference,
        "relative_cross_pair_count_difference": pair_count_difference,
        "pair_histogram_L1_fraction": hist_l1_fraction,
        "minimum_distance_difference_A": abs(a["minimum_Ag_W_distance_A"] - b["minimum_Ag_W_distance_A"]),
        "interpretation": "localized EAM per-atom energy is a diagnostic partition; the total cell energy remains the reported thermodynamic value.",
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lammps", default="/workspace/.venvs/ag-qsgs/bin/lmp")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if not SOURCE_REPORT.is_file():
        raise FileNotFoundError(SOURCE_REPORT)
    executable = str(args.lammps)
    report = {
        "title": "Independent audit of both periodic interfaces in the relaxed EFS screening cells",
        "potential_sha256": __import__("hashlib").sha256(POTENTIAL.read_bytes()).hexdigest(),
        "method": {
            "cells": "periodic x/y/z with two Ag/W boundaries per cell; bulk-like central cores fixed; same EFS minimized structures and 2.8 A initial gap as the reported lowest-energy registry",
            "locality_window_A": 2.0 * WINDOW_HALF_WIDTH_A,
            "per_atom_energy_partition": "LAMMPS pe/atom; subtract same-strain element bulk energy and divide by the lateral area",
            "cross_pair_distance_bins_A": DISTANCE_BINS_A,
            "caveat": "Local pe/atom windows are diagnostic, not a unique thermodynamic partition. The pair of interfaces is still thermodynamically reported as an average unless their local environments match.",
        },
        "DFT_reference_conditions_from_He_2018": {
            "method": "CASTEP PBE; 360 eV cutoff; bulk 9x9x9 k mesh; surfaces/interfaces 9x9x1; energy convergence 5e-5 eV/atom; force 0.03 eV/A; displacement 0.01 A",
            "surface_slabs": "7 atomic layers reported converged; two free surfaces with 15 A vacuum",
            "interface_cell": "interface slab in vacuum (one internal interface), Table 4 Wsep and gamma definitions",
            "Ag100_W100_orientation": "[010]Ag parallel [011]W; reported mismatch about 3.44%",
            "Ag111_W110_orientation": "[-2,1,1]Ag parallel [1,-1,0]W; reported mismatch about 3.44%",
        },
        "orientations": {},
    }
    with tempfile.TemporaryDirectory(prefix="agw-interface-equivalence-") as temporary:
        work = Path(temporary)
        for name, case in CASES.items():
            report["orientations"][name] = run_minimized_case(
                executable, name, case["registry"], case["gap_A"], work
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
