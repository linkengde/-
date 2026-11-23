#!/usr/bin/env python3
"""Relax small periodic Ag/W interface cells with pinned bulk-like cores."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell
from ase.io import write


ROOT = Path(__file__).resolve().parents[1]
POTENTIAL = ROOT / "potential_files" / "derived" / "Dai_2006_EFS_AgW.eam.fs"
OUTPUT = ROOT / "validation" / "relaxed_interface_report.json"
CONVERSION_J_M2 = 16.021766208
GAPS_A = (2.0, 2.4, 2.8, 3.2, 3.6, 6.0)
VACUUM_A = 12.0

INTERFACES = {
    "Ag100_W100": {
        "ag_matrix": [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
        "w_matrix": [[0, 1, 1], [0, -1, 1], [1, 0, 0]],
        "ag_repeats": (2, 2, 4),
        "w_repeats": (2, 2, 7),
        "core_depth_A": 3.0,
        "registries": {
            "top": (0.0, 0.0),
            "bridge_x": (0.25, 0.0),
            "bridge_diagonal": (0.25, 0.125),
        },
        "dft": {"work_of_separation_J_m2": 2.80, "interface_energy_J_m2": 0.72},
    },
    "Ag111_W110": {
        "ag_matrix": [[-2, 1, 1], [0, -1, 1], [1, 1, 1]],
        "w_matrix": [[1, -1, 0], [0, 0, -1], [1, 1, 0]],
        "ag_repeats": (3, 5, 3),
        "w_repeats": (7, 9, 4),
        "core_depth_A": 3.5,
        "registries": {
            "origin": (0.0, 0.0),
            "shift_x": (0.25, 0.0),
            "shift_xy": (0.125, 0.125),
            "shift_mixed": (0.375, 0.25),
        },
        "dft": {"work_of_separation_J_m2": 2.45, "interface_energy_J_m2": 0.86},
    },
}


def orient_cell_axes(atoms: Atoms) -> Atoms:
    rotation = (atoms.cell.array / atoms.cell.lengths()[:, None]).T
    atoms.positions = atoms.positions @ rotation
    atoms.set_cell(atoms.cell.array @ rotation)
    return atoms


def make_phase_slab(
    element: str,
    matrix: list[list[int]],
    repeats: tuple[int, int, int],
    target_xy: tuple[float, float],
) -> Atoms:
    structure = "fcc" if element == "Ag" else "bcc"
    lattice = 4.09 if element == "Ag" else 3.16
    primitive = orient_cell_axes(make_supercell(bulk(element, structure, a=lattice, cubic=True), matrix))
    slab = primitive.repeat(repeats)
    positions = slab.positions.copy()
    original_lengths = slab.cell.lengths()
    positions[:, 0] *= target_xy[0] / original_lengths[0]
    positions[:, 1] *= target_xy[1] / original_lengths[1]
    cell = slab.cell.array.copy()
    cell[0, 0] = target_xy[0]
    cell[1, 1] = target_xy[1]
    slab.positions = positions
    slab.set_cell(cell)
    return slab


def build_orientation(name: str) -> dict:
    spec = INTERFACES[name]
    ag_lattice = 4.09
    w_lattice = 3.16
    ag_primitive = orient_cell_axes(
        make_supercell(bulk("Ag", "fcc", a=ag_lattice, cubic=True), spec["ag_matrix"])
    )
    w_primitive = orient_cell_axes(
        make_supercell(bulk("W", "bcc", a=w_lattice, cubic=True), spec["w_matrix"])
    )
    ag_repeats = spec["ag_repeats"]
    w_repeats = spec["w_repeats"]
    ag_natural = ag_primitive.cell.lengths()[:2] * np.array(ag_repeats[:2])
    w_natural = w_primitive.cell.lengths()[:2] * np.array(w_repeats[:2])
    common_xy = np.sqrt(ag_natural * w_natural)
    ag = make_phase_slab("Ag", spec["ag_matrix"], ag_repeats, common_xy)
    w = make_phase_slab("W", spec["w_matrix"], w_repeats, common_xy)
    ag_strain = common_xy / ag_natural - 1.0
    w_strain = common_xy / w_natural - 1.0
    area = float(np.prod(common_xy))
    return {
        "ag": ag,
        "w": w,
        "area_A2": area,
        "lateral_size_A": common_xy,
        "ag_lateral_strain": ag_strain,
        "w_lateral_strain": w_strain,
    }


def core_region(bounds: tuple[float, float], name: str) -> str:
    low, high = bounds
    if low >= high:
        raise ValueError(f"Not enough layers to define a pinned {name} core")
    return (
        f"region {name} block INF INF INF INF "
        f"{low:.10f} {high:.10f} units box\n"
        f"group {name} region {name}\n"
    )


def build_interface_atoms(
    geometry: dict,
    gap: float,
    registry_fraction: tuple[float, float],
    core_depth: float,
) -> tuple[Atoms, dict]:
    ag = geometry["ag"]
    w = geometry["w"]
    lateral = geometry["lateral_size_A"]
    ag_positions = ag.positions.copy()
    w_positions = w.positions.copy()
    w_positions[:, :2] += lateral * np.asarray(registry_fraction)
    w_positions[:, :2] %= lateral

    ag_min, ag_max = float(ag_positions[:, 2].min()), float(ag_positions[:, 2].max())
    w_min, w_max = float(w_positions[:, 2].min()), float(w_positions[:, 2].max())
    w_positions[:, 2] += ag_max + gap - w_min
    w_low = ag_max + gap
    w_high = ag_max + gap + w_max - w_min
    box_z = w_high + gap - ag_min
    ag_positions[:, 2] -= ag_min

    positions = np.vstack((ag_positions, w_positions))
    symbols = ["Ag"] * len(ag) + ["W"] * len(w)
    atoms = Atoms(symbols=symbols, positions=positions, cell=[*lateral, box_z], pbc=True)
    atom_symbols = np.asarray(symbols)
    fixed_regions = {
        "ag_core": (core_depth, ag_max - ag_min - core_depth),
        "w_core": (w_low + core_depth, w_high - core_depth),
    }
    counts = {
        "Ag": int(np.count_nonzero(atom_symbols == "Ag")),
        "W": int(np.count_nonzero(atom_symbols == "W")),
        "fixed_core_atoms": 0,
    }
    for element, bounds in fixed_regions.items():
        low, high = bounds
        if low < high:
            if element == "ag_core":
                core_count = np.count_nonzero((ag_positions[:, 2] >= low) & (ag_positions[:, 2] <= high))
            else:
                core_count = np.count_nonzero((w_positions[:, 2] >= low) & (w_positions[:, 2] <= high))
            counts["fixed_core_atoms"] += int(core_count)
    if counts["fixed_core_atoms"] == 0:
        raise ValueError("No atoms fall in the pinned bulk-like core regions")
    return atoms, {"core_regions": fixed_regions, "counts": counts, "z_bounds": [0.0, box_z]}


def run_lammps(
    lammps: str,
    atoms: Atoms,
    work: Path,
    run_name: str,
    fixed_regions: dict[str, tuple[float, float]],
) -> dict:
    directory = work / run_name
    directory.mkdir()
    species = list(dict.fromkeys(atoms.get_chemical_symbols()))
    write(
        directory / "atoms.data",
        atoms,
        format="lammps-data",
        atom_style="atomic",
        masses=True,
        specorder=species,
    )
    regions = "".join(core_region(bounds, name) for name, bounds in fixed_regions.items())
    group_names = list(fixed_regions)
    if not group_names:
        raise ValueError(f"No fixed bulk-like core was defined for {run_name}")
    if len(group_names) == 1:
        core_group = f"group core region {group_names[0]}"
    else:
        core_group = f"group core union {' '.join(group_names)}"
    element_map = " ".join(species)
    (directory / "in.lmp").write_text(
        f"""units metal
atom_style atomic
boundary p p p
read_data atoms.data
pair_style eam/fs
pair_coeff * * {POTENTIAL} {element_map}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
{regions}{core_group}
fix hold core setforce 0.0 0.0 0.0
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo 100
thermo_style custom step atoms pe c_fmax
run 0
variable einitial equal pe
print "INITIAL energy=$(v_einitial:%.14g)"
min_style cg
minimize 0.0 1e-5 8000 80000
variable etotal equal pe
print "RESULT energy=$(v_etotal:%.14g) max_mobile_force=$(c_fmax:%.12g)"
"""
    )
    result = subprocess.run([lammps, "-in", "in.lmp", "-log", "none"], cwd=directory, capture_output=True, text=True)
    output = result.stdout + result.stderr
    if result.returncode:
        raise RuntimeError(f"LAMMPS failed for {run_name}:\n{output[-4000:]}")
    match = re.search(r"RESULT energy=([-+0-9.eE]+) max_mobile_force=([-+0-9.eE]+)", output)
    if match is None:
        raise RuntimeError(f"LAMMPS did not report minimized results for {run_name}:\n{output[-2500:]}")
    initial_match = re.search(r"INITIAL energy=([-+0-9.eE]+)", output)
    stopping_matches = re.findall(r"Stopping criterion = (.+)", output)
    iteration_matches = re.findall(r"Iterations, force evaluations =\s*(\d+)\s+\d+", output)
    return {
        "initial_energy_eV": float(initial_match.group(1)) if initial_match else None,
        "energy_eV": float(match.group(1)),
        "max_mobile_force_eV_A": float(match.group(2)),
        "minimization_iterations": int(iteration_matches[-1]) if iteration_matches else None,
        "stopping_criterion": stopping_matches[-1].strip() if stopping_matches else None,
    }


def surface_references(
    lammps: str,
    name: str,
    geometry: dict,
    core_depth: float,
    work: Path,
) -> dict:
    results = {}
    area = geometry["area_A2"]
    for element in ("Ag", "W"):
        slab = geometry[element.lower()]
        positions = slab.positions.copy()
        zmin, zmax = float(positions[:, 2].min()), float(positions[:, 2].max())
        slab_box = slab.cell.array.copy()
        slab_box[2, 2] += VACUUM_A
        surface_atoms = Atoms(
            symbols=[element] * len(slab),
            positions=positions,
            cell=slab_box,
            pbc=True,
        )
        if element == "Ag":
            regions = {"ag_core": (zmin + core_depth, zmax - core_depth)}
        else:
            regions = {"w_core": (zmin + core_depth, zmax - core_depth)}
        surface_result = run_lammps(lammps, surface_atoms, work, f"{name}_{element}_surface", regions)
        surface_energy = surface_result["energy_eV"]

        # The periodic reference keeps the interface's in-plane strain and relaxes its normal repeat length.
        primitive_matrix = INTERFACES[name][f"{element.lower()}_matrix"]
        repeats = INTERFACES[name][f"{element.lower()}_repeats"]
        structure = "fcc" if element == "Ag" else "bcc"
        lattice = 4.09 if element == "Ag" else 3.16
        primitive = orient_cell_axes(make_supercell(bulk(element, structure, a=lattice, cubic=True), primitive_matrix))
        bulk_atoms = primitive.repeat(repeats)
        bulk_positions = bulk_atoms.positions.copy()
        original_lengths = bulk_atoms.cell.lengths()
        target_xy = geometry["lateral_size_A"]
        bulk_positions[:, 0] *= target_xy[0] / original_lengths[0]
        bulk_positions[:, 1] *= target_xy[1] / original_lengths[1]
        bulk_cell = bulk_atoms.cell.array.copy()
        bulk_cell[0, 0] = target_xy[0]
        bulk_cell[1, 1] = target_xy[1]
        bulk_atoms.positions = bulk_positions
        bulk_atoms.set_cell(bulk_cell)
        bulk_directory = work / f"{name}_{element}_bulk"
        bulk_directory.mkdir()
        write(
            bulk_directory / "atoms.data",
            bulk_atoms,
            format="lammps-data",
            atom_style="atomic",
            masses=True,
            specorder=[element],
        )
        (bulk_directory / "in.lmp").write_text(
            f"""units metal
atom_style atomic
boundary p p p
read_data atoms.data
pair_style eam/fs
pair_coeff * * {POTENTIAL} {element}
neighbor 2.0 bin
fix relax all box/relax z 0.0 vmax 0.001
min_style cg
minimize 1e-10 1e-10 2000 20000
variable epa equal pe/count(all)
print "RESULT bulk_epa=$(v_epa:%.14g)"
"""
        )
        result = subprocess.run(
            [lammps, "-in", "in.lmp", "-log", "none"],
            cwd=bulk_directory,
            capture_output=True,
            text=True,
        )
        output = result.stdout + result.stderr
        if result.returncode:
            raise RuntimeError(f"LAMMPS failed for {name} {element} bulk: {output[-3000:]}")
        match = re.search(r"RESULT bulk_epa=([-+0-9.eE]+)", output)
        if match is None:
            raise RuntimeError(f"LAMMPS did not report {name} {element} bulk energy")
        bulk_epa = float(match.group(1))
        surface_gamma = (surface_energy - len(slab) * bulk_epa) / (2.0 * area) * CONVERSION_J_M2
        results[element] = {
            "bulk_energy_eV_atom": bulk_epa,
            "surface_slab_energy_eV": surface_energy,
            "surface_energy_J_m2": surface_gamma,
            "surface_atoms": len(slab),
            "surface_mobile_max_force_eV_A": surface_result["max_mobile_force_eV_A"],
        }
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lammps", default=os.environ.get("LAMMPS_EXECUTABLE", "/workspace/.venvs/ag-qsgs/bin/lmp"))
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--orientations", nargs="+", choices=tuple(INTERFACES), default=list(INTERFACES))
    args = parser.parse_args()
    lammps = shutil.which(args.lammps) or (args.lammps if Path(args.lammps).is_file() else None)
    if not lammps:
        raise FileNotFoundError(f"LAMMPS executable not found: {args.lammps}")
    if not POTENTIAL.is_file():
        raise FileNotFoundError(f"Generate the potential first: {POTENTIAL}")

    report = {
        "title": "Relaxed Dai EFS Ag/W interface scans",
        "potential_sha256": __import__("hashlib").sha256(POTENTIAL.read_bytes()).hexdigest(),
        "lammps_version": subprocess.run([lammps, "-help"], capture_output=True, text=True, check=True).stdout.splitlines()[1].strip(),
        "protocol": {
            "method": "LAMMPS energy minimization at fixed in-plane cell and initial vertical slab offset",
            "core_constraint": "bulk-like central slab atoms fixed with zero force; interfacial atoms relaxed",
            "interfaces_per_periodic_cell": 2,
            "gaps_A": list(GAPS_A),
            "separated_reference_gap_A": max(GAPS_A),
            "surface_vacuum_A": VACUUM_A,
            "minimizer": "LAMMPS conjugate-gradient, energy tolerance 0, force tolerance 1e-5 eV/A, maximum 8000 iterations and 80000 force evaluations",
            "work_of_separation": "(E_Ag_surface + E_W_surface - E_interface)/(2*area)",
            "interface_energy": "(E_interface - N_Ag*e_Ag_bulk - N_W*e_W_bulk)/(2*area)",
        },
        "interpretation_limitations": [
            "This is a zero-temperature EFS minimization benchmark, not a DFT-quality convergence study.",
            "Lateral cell dimensions are fixed to the selected commensurate match; the reported phase strains are part of the comparison.",
            "Bulk-like slab cores are held fixed while the remaining atoms relax; the cell height is fixed during each gap point.",
            "Each periodic supercell contains two Ag/W interfaces. No production or finite-temperature MD is performed.",
        ],
        "orientations": {},
    }

    with tempfile.TemporaryDirectory(prefix="ag-w-relaxed-interface-") as temporary:
        work = Path(temporary)
        for name in args.orientations:
            spec = INTERFACES[name]
            geometry = build_orientation(name)
            area = geometry["area_A2"]
            surfaces = surface_references(lammps, name, geometry, spec["core_depth_A"], work)
            orientation_result = {
                "dft_reference": spec["dft"],
                "area_A2": area,
                "lateral_size_A": geometry["lateral_size_A"].tolist(),
                "Ag_lateral_strain": geometry["ag_lateral_strain"].tolist(),
                "W_lateral_strain": geometry["w_lateral_strain"].tolist(),
                "surface_references": surfaces,
                "registries": {},
            }
            for registry, fraction in spec["registries"].items():
                scan = []
                for gap in GAPS_A:
                    atoms, build_info = build_interface_atoms(
                        geometry,
                        gap,
                        fraction,
                        spec["core_depth_A"],
                    )
                    result = run_lammps(
                        lammps,
                        atoms,
                        work,
                        f"{name}_{registry}_gap_{gap:.1f}",
                        build_info["core_regions"],
                    )
                    result.update(
                        {
                            "initial_vertical_gap_A": gap,
                            "initial_registry_fraction": list(fraction),
                            **build_info["counts"],
                            "box_z_A": float(atoms.cell.lengths()[2]),
                        }
                    )
                    scan.append(result)

                reference = next(item for item in scan if item["initial_vertical_gap_A"] == max(GAPS_A))
                minimum = min(scan, key=lambda item: item["energy_eV"])
                binding_per_interface = (minimum["energy_eV"] - reference["energy_eV"]) / (2.0 * area)
                wsep = -binding_per_interface * CONVERSION_J_M2
                ag = surfaces["Ag"]
                w = surfaces["W"]
                interface_energy = (
                    minimum["energy_eV"]
                    - ag["surface_atoms"] * ag["bulk_energy_eV_atom"]
                    - w["surface_atoms"] * w["bulk_energy_eV_atom"]
                ) / (2.0 * area) * CONVERSION_J_M2
                surface_energy_sum = ag["surface_energy_J_m2"] + w["surface_energy_J_m2"]
                scan_min = {
                    "minimum_energy_initial_gap_A": minimum["initial_vertical_gap_A"],
                    "minimum_energy_eV": minimum["energy_eV"],
                    "separated_reference_energy_eV": reference["energy_eV"],
                    "separated_reference_minus_surface_slabs_eV": reference["energy_eV"]
                    - ag["surface_slab_energy_eV"]
                    - w["surface_slab_energy_eV"],
                    "work_of_separation_J_m2": wsep,
                    "interface_energy_J_m2": interface_energy,
                    "surface_energy_sum_J_m2": surface_energy_sum,
                    "work_of_separation_closure_error_J_m2": surface_energy_sum - interface_energy - wsep,
                    "wsep_error_vs_dft_J_m2": wsep - spec["dft"]["work_of_separation_J_m2"],
                    "interface_energy_error_vs_dft_J_m2": interface_energy - spec["dft"]["interface_energy_J_m2"],
                    "max_mobile_force_at_minimum_eV_A": minimum["max_mobile_force_eV_A"],
                    "minimization_iterations_at_minimum": minimum["minimization_iterations"],
                }
                orientation_result["registries"][registry] = {"scan": scan, "summary": scan_min}
            report["orientations"][name] = orientation_result

    report["status"] = "RELAXED_INTERFACE_SCANS_COMPLETE_PRELIMINARY"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
