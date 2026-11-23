#!/usr/bin/env python3
"""Build and screen low-strain, one-interface Ag/W slabs with the Dai EFS potential."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell
from ase.neighborlist import neighbor_list
from ase.io import write


ROOT = Path(__file__).resolve().parents[3]
POTENTIAL = ROOT / "potential_files/derived/Dai_2006_EFS_AgW.eam.fs"
MATCH = Path(__file__).resolve().parent / "cell_match_report.json"
OUT = Path(__file__).resolve().parent
CONVERSION_J_M2 = 16.021766208
N_LAYERS = 7
N_FIXED_LAYERS = 3
VACUUM_A = 15.0
GAPS_A = tuple([round(2.0 + 0.2 * i, 1) for i in range(9)] + [6.0])
FORCE_TOL_EVA = 0.01
DFT_FORCE_LIMIT_EVA = 0.03
SAME_PHASE_COLLAPSE_FRACTION = 0.85

ORIENTATION = {
    "Ag100_W100": {
        "Ag_matrix": [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
        "W_matrix": [[0, 1, 1], [0, -1, 1], [1, 0, 0]],
        "registries": {
            "top": (0.0, 0.0),
            "bridge_x": (0.5, 0.0),
            "bcc_type": (0.5, 0.5),
        },
        "dft": {"work_of_separation_J_m2": 2.80, "interface_energy_J_m2": 0.72},
    },
    "Ag111_W110": {
        "Ag_matrix": [[-2, 1, 1], [0, -1, 1], [1, 1, 1]],
        "W_matrix": [[1, -1, 0], [0, 0, -1], [1, 1, 0]],
        "registries": {
            "shift_00": (0.0, 0.0),
            "shift_half_x": (0.5, 0.0),
            "shift_half_y": (0.0, 0.5),
            "shift_half_xy": (0.5, 0.5),
        },
        "dft": {"work_of_separation_J_m2": 2.45, "interface_energy_J_m2": 0.86},
    },
}


def orient_cell_axes(atoms: Atoms) -> Atoms:
    rotation = (atoms.cell.array / atoms.cell.lengths()[:, None]).T
    atoms.positions = atoms.positions @ rotation
    atoms.set_cell(atoms.cell.array @ rotation)
    return atoms


def cluster_levels(values: np.ndarray, tolerance: float = 1.0e-5) -> list[float]:
    levels: list[list[float]] = []
    for value in sorted(float(v) for v in values):
        if not levels or abs(value - float(np.mean(levels[-1]))) > tolerance:
            levels.append([value])
        else:
            levels[-1].append(value)
    return [float(np.mean(level)) for level in levels]


def oriented_bulk(element: str, matrix: list[list[int]], lattice: float) -> Atoms:
    crystal = "fcc" if element == "Ag" else "bcc"
    return orient_cell_axes(make_supercell(bulk(element, crystal, a=lattice, cubic=True), matrix))


def make_phase_slab(
    element: str,
    matrix: list[list[int]],
    repeats_xy: list[int],
    common_xy: list[float],
    lattice: float,
    termination: int = 0,
) -> tuple[Atoms, dict]:
    base = oriented_bulk(element, matrix, lattice)
    planes_per_cell = len(cluster_levels(base.positions[:, 2]))
    z_repeat = math.ceil((N_LAYERS + termination) / planes_per_cell) + 1
    stack = base.repeat((int(repeats_xy[0]), int(repeats_xy[1]), z_repeat))
    positions = stack.positions.copy()
    natural_x = float(base.cell.lengths()[0] * repeats_xy[0])
    natural_y = float(base.cell.lengths()[1] * repeats_xy[1])
    positions[:, 0] *= common_xy[0] / natural_x
    positions[:, 1] *= common_xy[1] / natural_y
    levels = cluster_levels(positions[:, 2])
    selected_levels = levels[termination : termination + N_LAYERS]
    if len(selected_levels) != N_LAYERS:
        raise ValueError(f"Could not select {N_LAYERS} planes for {element}; got {len(selected_levels)}")
    layer_id = np.full(len(stack), -1, dtype=int)
    for index, level in enumerate(selected_levels):
        layer_id[np.abs(positions[:, 2] - level) < 1.0e-5] = index
    mask = layer_id >= 0
    if not np.all(mask):
        # Atoms outside the seven selected planes are discarded by construction, not by relaxation.
        pass
    selected = positions[mask].copy()
    selected_layer_id = layer_id[mask]
    z_origin = float(selected_levels[0])
    selected[:, 2] -= z_origin
    thickness = float(selected[:, 2].max())
    cell_z = thickness + VACUUM_A
    atoms = Atoms(
        symbols=[element] * len(selected),
        positions=selected,
        cell=np.diag([common_xy[0], common_xy[1], cell_z]),
        pbc=(True, True, True),
    )
    atoms.new_array("layer_id", selected_layer_id)
    core_low = 0.5 * (selected_levels[1] + selected_levels[2]) - z_origin
    core_high = 0.5 * (selected_levels[4] + selected_levels[5]) - z_origin
    metadata = {
        "element": element,
        "termination_index": termination,
        "planes_per_oriented_repeat": planes_per_cell,
        "layer_z_A": (np.asarray(selected_levels) - z_origin).tolist(),
        "core_region_z_A": [float(core_low), float(core_high)],
        "material_thickness_A": thickness,
        "atom_count": len(atoms),
        "natural_lateral_lengths_A": [natural_x, natural_y],
        "common_lateral_lengths_A": list(common_xy),
        "applied_strain": [common_xy[0] / natural_x - 1.0, common_xy[1] / natural_y - 1.0],
    }
    return atoms, metadata


def build_interface(
    ag: Atoms,
    ag_meta: dict,
    w: Atoms,
    w_meta: dict,
    common_xy: list[float],
    translation_fraction: tuple[float, float],
    gap_A: float,
) -> tuple[Atoms, dict]:
    ag_pos = ag.positions.copy()
    w_pos = w.positions.copy()
    shift_xy = np.asarray(translation_fraction) * np.asarray(common_xy)
    w_pos[:, 0] = np.mod(w_pos[:, 0] + shift_xy[0], common_xy[0])
    w_pos[:, 1] = np.mod(w_pos[:, 1] + shift_xy[1], common_xy[1])
    ag_thick = float(ag_pos[:, 2].max())
    w_offset = ag_thick + gap_A
    w_pos[:, 2] += w_offset
    material_top = w_offset + float(w_meta["material_thickness_A"])
    positions = np.vstack((ag_pos, w_pos))
    symbols = ["Ag"] * len(ag) + ["W"] * len(w)
    atoms = Atoms(
        symbols=symbols,
        positions=positions,
        cell=np.diag([common_xy[0], common_xy[1], material_top + VACUUM_A]),
        pbc=(True, True, True),
    )
    ag_core = tuple(float(x) for x in ag_meta["core_region_z_A"])
    w_core = tuple(float(x) + w_offset for x in w_meta["core_region_z_A"])
    metadata = {
        "Ag_atoms": len(ag),
        "W_atoms": len(w),
        "atom_count": len(atoms),
        "Ag_core_region_z_A": list(ag_core),
        "W_core_region_z_A": list(w_core),
        "interface_gap_A": float(gap_A),
        "translation_fraction_of_common_cell": list(translation_fraction),
        "translation_A": shift_xy.tolist(),
        "fixed_core_atoms": int(np.count_nonzero((ag.arrays["layer_id"] >= 2) & (ag.arrays["layer_id"] <= 4)))
        + int(np.count_nonzero((w.arrays["layer_id"] >= 2) & (w.arrays["layer_id"] <= 4))),
    }
    return atoms, metadata


def lammps_groups(core_regions: dict[str, tuple[float, float]]) -> tuple[str, str]:
    lines = []
    names = []
    for name, (lo, hi) in core_regions.items():
        lines.extend(
            [
                f"region {name}_region block INF INF INF INF {lo:.12f} {hi:.12f} units box",
                f"group {name} region {name}_region",
            ]
        )
        names.append(name)
    if len(names) == 1:
        core = f"group core union {names[0]}"
    else:
        core = f"group core union {' '.join(names)}"
    lines.extend([core, "group mobile subtract all core"])
    return "\n".join(lines), "\n".join(names)


def run_lammps(
    lammps: str,
    atoms: Atoms,
    directory: Path,
    species_order: list[str],
    core_regions: dict[str, tuple[float, float]] | None = None,
    relax_bulk_z: bool = False,
) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    data_file = directory / "atoms.data"
    write(data_file, atoms, format="lammps-data", atom_style="atomic", masses=True, specorder=species_order)
    groups = ""
    if core_regions:
        groups, _ = lammps_groups(core_regions)
    potential_map = " ".join(species_order)
    if relax_bulk_z:
        relaxation = "fix relax all box/relax z 0.0 vmax 0.001\n"
        group_name = "all"
    else:
        relaxation = ""
        group_name = "mobile"
    input_text = f"""units metal
atom_style atomic
boundary p p p
read_data atoms.data
pair_style eam/fs
pair_coeff * * {POTENTIAL} {potential_map}
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
{groups}
{relaxation}variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute mobile_fmax {group_name} reduce max v_fmag
thermo 100
thermo_style custom step atoms pe c_mobile_fmax fmax
thermo_modify lost error
min_style cg
min_modify dmax 0.01 line quadratic
minimize 0.0 {FORCE_TOL_EVA:.8g} 10000 100000
run 0
variable initial_energy equal pe
print "INITIAL energy=$(v_initial_energy:%.14g)"
variable total_energy equal pe
print "RESULT energy=$(v_total_energy:%.14g) max_mobile_force=$(c_mobile_fmax:%.12g)"
write_data minimized.data
write_dump all custom minimized.dump id type x y z modify sort id
"""
    (directory / "in.lmp").write_text(input_text)
    proc = subprocess.run(
        [lammps, "-in", "in.lmp", "-log", "log.lammps"],
        cwd=directory,
        capture_output=True,
        text=True,
    )
    output = proc.stdout + proc.stderr
    if proc.returncode:
        return {"run_ok": False, "error": output[-5000:], "work_directory": str(directory)}
    result_match = re.search(r"RESULT energy=([-+0-9.eE]+) max_mobile_force=([-+0-9.eE]+)", output)
    if result_match is None:
        return {"run_ok": False, "error": f"LAMMPS result line not parsed\n{output[-3000:]}", "work_directory": str(directory)}
    stops = re.findall(r"Stopping criterion = (.+)", output)
    iterations = re.findall(r"Iterations, force evaluations =\s*(\d+)\s+\d+", output)
    initial_match = re.search(r"INITIAL energy=([-+0-9.eE]+)", output)
    return {
        "run_ok": True,
        "initial_energy_eV": float(initial_match.group(1)) if initial_match else None,
        "energy_eV": float(result_match.group(1)),
        "max_mobile_force_eV_A": float(result_match.group(2)),
        "stopping_criterion": stops[-1].strip() if stops else None,
        "iterations": int(iterations[-1]) if iterations else None,
        "atom_count": len(atoms),
        "minimized_data": str(directory / "minimized.data"),
        "minimized_dump": str(directory / "minimized.dump"),
        "work_directory": str(directory),
    }


def report_copy_without_temporary_paths(value):
    if isinstance(value, dict):
        return {
            key: report_copy_without_temporary_paths(item)
            for key, item in value.items()
            if key not in {"minimized_data", "minimized_dump", "work_directory"}
        }
    if isinstance(value, list):
        return [report_copy_without_temporary_paths(item) for item in value]
    return value


def save_report(path: Path, report: dict) -> None:
    path.write_text(json.dumps(report_copy_without_temporary_paths(report), indent=2) + "\n")


def read_minimized_dump(path: Path, template: Atoms, species_order: list[str]) -> Atoms:
    lines = path.read_text().splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS"))
    fields = lines[header].split()[2:]
    indexes = {name: fields.index(name) for name in ("id", "type", "x", "y", "z")}
    rows = [line.split() for line in lines[header + 1 : header + 1 + len(template)]]
    rows.sort(key=lambda row: int(row[indexes["id"]]))
    types = np.array([int(row[indexes["type"]]) for row in rows], dtype=int)
    symbols_by_type = {i + 1: symbol for i, symbol in enumerate(species_order)}
    symbols = [symbols_by_type[int(t)] for t in types]
    positions = np.array([[float(row[indexes[k]]) for k in ("x", "y", "z")] for row in rows])
    return Atoms(symbols=symbols, positions=positions, cell=template.cell.array, pbc=template.pbc)


def distance_summary(atoms: Atoms) -> dict:
    i, j, d = neighbor_list("ijd", atoms, 4.76, self_interaction=False)
    sym = np.asarray(atoms.get_chemical_symbols())
    values = {"Ag_Ag": [], "W_W": [], "Ag_W": []}
    for ii, jj, dd in zip(i, j, d):
        if ii >= jj:
            continue
        pair = tuple(sorted((sym[ii], sym[jj])))
        key = "Ag_Ag" if pair == ("Ag", "Ag") else "W_W" if pair == ("W", "W") else "Ag_W"
        values[key].append(float(dd))
    return {key: min(items) if items else None for key, items in values.items()}


def compress_file(path: Path, target: Path) -> None:
    with path.open("rb") as src, target.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)


def build_orientation_geometry(name: str, match: dict, termination_ag: int = 0, termination_w: int = 0):
    spec = ORIENTATION[name]
    cell = match["orientations"][name]["selected_cell"]
    common_xy = [float(v) for v in cell["common_lengths_A"]]
    lattices = match["equilibrium_lattice_A"]
    ag, ag_meta = make_phase_slab(
        "Ag", spec["Ag_matrix"], cell["Ag_repeats_xy"], common_xy, lattices["Ag"], termination_ag
    )
    w, w_meta = make_phase_slab(
        "W", spec["W_matrix"], cell["W_repeats_xy"], common_xy, lattices["W"], termination_w
    )
    return ag, ag_meta, w, w_meta, common_xy


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lammps", default=os.environ.get("LAMMPS_EXECUTABLE", "/workspace/.venvs/ag-qsgs/bin/lmp"))
    parser.add_argument("--orientations", nargs="+", choices=tuple(ORIENTATION), default=list(ORIENTATION))
    parser.add_argument("--dry-run", action="store_true", help="build slabs and report counts without LAMMPS")
    parser.add_argument("--termination-ag", type=int, default=0)
    parser.add_argument("--termination-w", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUT / "local_interface_v3_report.json")
    args = parser.parse_args()
    lammps = shutil.which(args.lammps) or (args.lammps if Path(args.lammps).is_file() else None)
    if not lammps and not args.dry_run:
        raise FileNotFoundError(f"LAMMPS executable not found: {args.lammps}")
    if not POTENTIAL.is_file():
        raise FileNotFoundError(f"Dai EFS potential not found: {POTENTIAL}")
    match = json.loads(MATCH.read_text())
    version = None
    if lammps and not args.dry_run:
        help_output = subprocess.run([lammps, "-help"], capture_output=True, text=True, check=True).stdout
        version = help_output.splitlines()[1].strip()

    report = {
        "title": "Ag-W v3 one-interface low-strain local slab screening",
        "status": "GEOMETRY_PREFLIGHT_ONLY" if args.dry_run else "RUNNING",
        "potential_sha256": hashlib.sha256(POTENTIAL.read_bytes()).hexdigest(),
        "lammps_version": version,
        "frozen_protocol": {
            "layers_per_phase": N_LAYERS,
            "fixed_core_layers_per_phase": N_FIXED_LAYERS,
            "vacuum_A": VACUUM_A,
            "gaps_A": list(GAPS_A),
            "force_minimizer_tolerance_eV_A": FORCE_TOL_EVA,
            "pass_force_limit_eV_A": DFT_FORCE_LIMIT_EVA,
            "same_phase_collapse_fraction": SAME_PHASE_COLLAPSE_FRACTION,
            "one_interface_per_slab": True,
        },
        "dft_comparison_note": "He et al. report Table 4 values and approximately 3.44% mismatch but do not publish the full interface vectors/strain partition; this is a screening comparison, not a direct DFT reproduction.",
        "orientations": {},
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    structure_dir = OUT / "relaxed_structures"
    if not args.dry_run:
        structure_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="agw-v3-local-interface-") as temp_name:
        work = Path(temp_name)
        for name in args.orientations:
            spec = ORIENTATION[name]
            ag, ag_meta, w, w_meta, common_xy = build_orientation_geometry(
                name, match, args.termination_ag, args.termination_w
            )
            item = {
                "selected_match": match["orientations"][name]["selected_cell"],
                "termination_indices": {"Ag": args.termination_ag, "W": args.termination_w},
                "phase_slab_geometry": {"Ag": ag_meta, "W": w_meta},
                "dft_reference": spec["dft"],
                "atom_counts": {"Ag": len(ag), "W": len(w), "total": len(ag) + len(w)},
                "registries": {},
            }
            if args.dry_run:
                report["orientations"][name] = item
                continue

            area = common_xy[0] * common_xy[1]
            refs = {}
            for element, slab, metadata in (("Ag", ag, ag_meta), ("W", w, w_meta)):
                core = {f"{element.lower()}_core": tuple(metadata["core_region_z_A"])}
                surface_result = run_lammps(
                    lammps, slab, work / f"{name}_{element}_surface", [element], core_regions=core
                )
                if not surface_result["run_ok"]:
                    raise RuntimeError(surface_result["error"])
                bulk = oriented_bulk(element, spec[f"{element}_matrix"], match["equilibrium_lattice_A"][element])
                repeats_xy = item["selected_match"][f"{element}_repeats_xy"]
                bulk = bulk.repeat((int(repeats_xy[0]), int(repeats_xy[1]), 1))
                bulk_pos = bulk.positions.copy()
                bulk_pos[:, 0] *= common_xy[0] / float(bulk.cell.lengths()[0])
                bulk_pos[:, 1] *= common_xy[1] / float(bulk.cell.lengths()[1])
                bulk_cell = np.diag([common_xy[0], common_xy[1], float(bulk.cell.lengths()[2])])
                bulk_atoms = Atoms(
                    symbols=[element] * len(bulk),
                    positions=bulk_pos,
                    cell=bulk_cell,
                    pbc=(True, True, True),
                )
                bulk_result = run_lammps(
                    lammps, bulk_atoms, work / f"{name}_{element}_bulk", [element], relax_bulk_z=True
                )
                if not bulk_result["run_ok"]:
                    raise RuntimeError(bulk_result["error"])
                refs[element] = {
                    "bulk_energy_eV_atom": bulk_result["energy_eV"] / len(bulk_atoms),
                    "bulk_atom_count": len(bulk_atoms),
                    "surface_energy_J_m2": (surface_result["energy_eV"] - len(slab) * bulk_result["energy_eV"] / len(bulk_atoms))
                    / (2.0 * area)
                    * CONVERSION_J_M2,
                    "surface_energy_slab_eV": surface_result["energy_eV"],
                    "surface_atom_count": len(slab),
                    "surface_minimize": surface_result,
                    "bulk_minimize": bulk_result,
                }
            item["references"] = refs
            registry_best_structures = []
            for registry, fraction in spec["registries"].items():
                scan = []
                for gap in GAPS_A:
                    atoms, metadata = build_interface(ag, ag_meta, w, w_meta, common_xy, fraction, gap)
                    case = work / f"{name}_{registry}_gap_{gap:.1f}"
                    result = run_lammps(lammps, atoms, case, ["Ag", "W"], core_regions={
                        "ag_core": tuple(metadata["Ag_core_region_z_A"]),
                        "w_core": tuple(metadata["W_core_region_z_A"]),
                    })
                    if result["run_ok"]:
                        relaxed = read_minimized_dump(Path(result["minimized_dump"]), atoms, ["Ag", "W"])
                        result["minimum_distances_A"] = distance_summary(relaxed)
                        result["retained_atom_count"] = len(relaxed)
                        result["relaxed_cell_z_A"] = float(relaxed.cell.lengths()[2])
                    result.update(metadata)
                    result["initial_gap_A"] = float(gap)
                    scan.append(result)
                    print(f"[{name}] {registry} gap={gap:.1f} -> {result.get('stopping_criterion', result.get('error', 'failed'))}", flush=True)
                    item["registries"][registry] = {"scan": scan}
                    report["orientations"][name] = item
                    report["status"] = "RUNNING"
                    save_report(args.output, report)

                valid = [r for r in scan if r.get("run_ok")]
                if not valid:
                    continue
                separated = next(r for r in valid if r["initial_gap_A"] == 6.0)
                minimum = min(valid, key=lambda r: r["energy_eV"])
                ag_ref = refs["Ag"]
                w_ref = refs["W"]
                natoms = {"Ag": len(ag), "W": len(w)}
                bulk_sum = natoms["Ag"] * ag_ref["bulk_energy_eV_atom"] + natoms["W"] * w_ref["bulk_energy_eV_atom"]
                gamma_surfaces = ag_ref["surface_energy_J_m2"] + w_ref["surface_energy_J_m2"]
                wsep = (ag_ref["surface_energy_slab_eV"] + w_ref["surface_energy_slab_eV"] - minimum["energy_eV"]) / area * CONVERSION_J_M2
                gamma_int = (minimum["energy_eV"] - bulk_sum) / area * CONVERSION_J_M2 - gamma_surfaces
                closure = gamma_surfaces - gamma_int - wsep
                point_distances = minimum.get("minimum_distances_A", {})
                ref_nn = {"Ag_Ag": 4.0900219316 / np.sqrt(2.0), "W_W": np.sqrt(3.0) * 3.1599967891 / 2.0}
                collapse = any(
                    point_distances.get(pair) is not None and point_distances[pair] < SAME_PHASE_COLLAPSE_FRACTION * ref_nn[pair]
                    for pair in ("Ag_Ag", "W_W")
                )
                summary = {
                    "minimum_energy_initial_gap_A": minimum["initial_gap_A"],
                    "minimum_energy_eV": minimum["energy_eV"],
                    "separated_reference_energy_eV": separated["energy_eV"],
                    "separated_energy_minus_independent_surface_slabs_eV": separated["energy_eV"]
                    - ag_ref["surface_energy_slab_eV"]
                    - w_ref["surface_energy_slab_eV"],
                    "work_of_separation_J_m2": wsep,
                    "interface_energy_J_m2": gamma_int,
                    "sum_surface_energy_J_m2": gamma_surfaces,
                    "thermodynamic_closure_error_J_m2": closure,
                    "work_of_separation_error_vs_DFT_J_m2": wsep - spec["dft"]["work_of_separation_J_m2"],
                    "interface_energy_error_vs_DFT_J_m2": gamma_int - spec["dft"]["interface_energy_J_m2"],
                    "max_mobile_force_eV_A": minimum["max_mobile_force_eV_A"],
                    "minimization_stopping_criterion": minimum["stopping_criterion"],
                    "minimization_iterations": minimum["iterations"],
                    "minimum_distances_A": point_distances,
                    "same_phase_collapse_flag": collapse,
                    "force_gate_pass": minimum["stopping_criterion"] == "force tolerance"
                    and minimum["max_mobile_force_eV_A"] <= DFT_FORCE_LIMIT_EVA,
                    "atom_count_preserved": minimum.get("retained_atom_count") == len(ag) + len(w),
                }
                registry_best_structures.append((registry, summary, minimum))
                item["registries"][registry]["summary"] = summary
                best_data = Path(minimum["minimized_data"])
                compressed = structure_dir / (
                    f"{name}_tAg{args.termination_ag}_tW{args.termination_w}_"
                    f"{registry}_best_gap_{minimum['initial_gap_A']:.1f}.data.gz"
                )
                compress_file(best_data, compressed)
                summary["relaxed_structure_file"] = str(compressed.relative_to(OUT))
                summary["relaxed_structure_sha256"] = hashlib.sha256(compressed.read_bytes()).hexdigest()
                item["registries"][registry]["summary"] = summary

            if registry_best_structures:
                best_registry, best_summary, _ = min(registry_best_structures, key=lambda row: row[1]["minimum_energy_eV"])
                dft = spec["dft"]
                dft_agree = (
                    abs(best_summary["work_of_separation_error_vs_DFT_J_m2"]) <= 0.30
                    and abs(best_summary["interface_energy_error_vs_DFT_J_m2"]) <= 0.30
                    and best_summary["interface_energy_J_m2"] > 0
                )
                item["best_registry"] = best_registry
                item["best_interface_summary"] = best_summary
                item["mechanical_gate_pass"] = best_summary["force_gate_pass"] and not best_summary["same_phase_collapse_flag"]
                item["DFT_screening_gate_pass"] = dft_agree
                item["orientation_gate_pass"] = item["mechanical_gate_pass"] and dft_agree and best_summary["atom_count_preserved"]
            report["orientations"][name] = item
            report["status"] = "RUNNING"
            save_report(args.output, report)

    if args.dry_run:
        report["status"] = "GEOMETRY_PREFLIGHT_COMPLETE_NO_LAMMPS_RUN"
    else:
        passed = all(report["orientations"].get(name, {}).get("orientation_gate_pass", False) for name in args.orientations)
        report["status"] = "LOCAL_INTERFACE_SCREEN_PASSED" if passed else "LOCAL_INTERFACE_SCREEN_REJECTED"
    save_report(args.output, report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
