#!/usr/bin/env python3
"""Run pure-crystal and Ag/W interface checks for the generated EFS table."""

from __future__ import annotations

import argparse
import hashlib
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
from ase.data import atomic_masses, atomic_numbers
from ase.io import write
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
POTENTIAL = ROOT / "potential_files" / "derived" / "Dai_2006_EFS_AgW.eam.fs"
OUTPUT = ROOT / "validation" / "dai_2006_efs_validation_report.json"
PARAMETERS = ROOT / "parameters" / "Dai_2006_EFS_AgW_parameters.json"
REFERENCE = {
    "Ag": {"structure": "fcc", "lattice_A": 4.09, "cohesive_eV_atom": 2.95},
    "W": {"structure": "bcc", "lattice_A": 3.16, "cohesive_eV_atom": 8.916},
}
ALLOY_REFERENCE = {
    "AgW_B2": {
        "symbols": ["Ag", "W"],
        "scaled_positions": [[0, 0, 0], [0.5, 0.5, 0.5]],
        "lattice_A": 3.22,
        "cohesive_this_work_eV_atom": 4.7990,
        "cohesive_vasp_eV_atom": 4.7212,
    },
    "Ag3W_L12": {
        "symbols": ["W", "Ag", "Ag", "Ag"],
        "scaled_positions": [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
        "lattice_A": 4.13,
        "cohesive_this_work_eV_atom": 3.4291,
        "cohesive_vasp_eV_atom": 3.4542,
    },
}
GAPS_A = (1.5, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.6, 6.0)
TABLE_PARAMETERS = {
    "Ag": {
        "A": 0.325514, "d_A": 4.41, "c_A": 4.76, "c0": 10.6812,
        "c1": -12.04517, "c2": 5.203072, "c3": -1.013304,
        "c4": 0.0742308, "B": -1.293394,
    },
    "W": {
        "A": 1.885948, "d_A": 4.41, "c_A": 3.25, "c0": 48.52796,
        "c1": -33.79621, "c2": 5.854334, "c3": -0.0098221,
        "c4": 0.033338, "B": 0,
    },
    "Ag_W": {
        "A": 0.6962, "d_A": 4.6, "c_A": 4.3, "c0": 25.618,
        "c1": -26.3407, "c2": 9.1414, "c3": -1.0693,
        "c4": 0, "B": 0,
    },
}


def lammps_result(command: list[str], working_directory: Path, result_pattern: str) -> tuple[re.Match, str]:
    result = subprocess.run(command, cwd=working_directory, check=False, capture_output=True, text=True)
    output = result.stdout + result.stderr
    if result.returncode:
        raise RuntimeError(f"LAMMPS failed with code {result.returncode}:\n{output[-4000:]}")
    match = re.search(result_pattern, output)
    if match is None:
        raise RuntimeError(f"LAMMPS did not print its expected result:\n{output[-3000:]}")
    return match, output


def run_pure_crystal(lammps: str, element: str, work: Path) -> dict:
    benchmark = REFERENCE[element]
    directory = work / f"pure_{element}"
    directory.mkdir()
    cell = bulk(element, benchmark["structure"], a=benchmark["lattice_A"], cubic=True).repeat((4, 4, 4))
    write(directory / "atoms.data", cell, format="lammps-data", atom_style="atomic", masses=True, specorder=[element])
    (directory / "in.lmp").write_text(
        f"""units metal
atom_style atomic
boundary p p p
read_data atoms.data
pair_style eam/fs
pair_coeff * * {POTENTIAL} {element}
neighbor 2.0 bin
thermo 10
thermo_style custom step atoms pe press lx ly lz
fix relax all box/relax iso 0.0 vmax 0.001
min_style cg
minimize 1e-12 1e-12 1000 10000
unfix relax
variable epa equal pe/count(all)
variable lat equal lx/4.0
print "RESULT element {element} epa=$(v_epa:%.12g) lattice=$(v_lat:%.12g)"
"""
    )
    match, _ = lammps_result(
        [lammps, "-in", "in.lmp", "-log", "none"],
        directory,
        rf"RESULT element {element} epa=([-+0-9.eE]+) lattice=([-+0-9.eE]+)",
    )
    cohesive = -float(match.group(1))
    lattice = float(match.group(2))
    return {
        "structure": benchmark["structure"],
        "cohesive_eV_atom": cohesive,
        "cohesive_reference_eV_atom": benchmark["cohesive_eV_atom"],
        "cohesive_absolute_error_eV_atom": abs(cohesive - benchmark["cohesive_eV_atom"]),
        "lattice_A": lattice,
        "lattice_reference_A": benchmark["lattice_A"],
        "lattice_absolute_error_A": abs(lattice - benchmark["lattice_A"]),
        "minimized_from_reported_lattice": True,
    }


def get_lammps_version(lammps: str) -> str:
    result = subprocess.run([lammps, "-help"], capture_output=True, text=True, check=True)
    output = result.stdout + result.stderr
    return next(line.strip() for line in output.splitlines() if "Atomic/Molecular Massively Parallel Simulator" in line)


def audit_source_parameters() -> dict:
    parameters = json.loads(PARAMETERS.read_text())
    if parameters["parameter_units"]["A"] != "eV Angstrom^-1":
        raise AssertionError("Dai Table 1, 4, and 6 list A in eV Angstrom^-1")
    groups = {
        "Ag": parameters["pure_element_parameters"]["Ag"],
        "W": parameters["pure_element_parameters"]["W"],
        "Ag_W": parameters["cross_parameters"]["Ag_W"],
    }
    for name, expected in TABLE_PARAMETERS.items():
        for field, value in expected.items():
            if groups[name][field] != value:
                raise AssertionError(f"{name}.{field} disagrees with the values checked in the Dai PDF")
    return {"source_parameters_sha256": hashlib.sha256(PARAMETERS.read_bytes()).hexdigest(), "table_values_match_audited_transcription": True}


def direct_efs_energy(atoms: Atoms) -> float:
    parameters = json.loads(PARAMETERS.read_text())
    symbols = np.asarray(atoms.get_chemical_symbols())
    positions = atoms.positions
    box = atoms.cell.lengths()
    maximum_cutoff = max(
        item["c_A"]
        for item in (
            *parameters["pure_element_parameters"].values(),
            *parameters["cross_parameters"].values(),
        )
    )
    pairs = cKDTree(positions, boxsize=box).query_pairs(maximum_cutoff, output_type="ndarray")
    density = np.zeros(len(atoms))
    pair_energy = 0.0
    for first, second in pairs:
        displacement = positions[second] - positions[first]
        displacement -= np.rint(displacement / box) * box
        radius = float(np.linalg.norm(displacement))
        first_symbol = symbols[first]
        second_symbol = symbols[second]
        pair = parameters["pure_element_parameters"][first_symbol] if first_symbol == second_symbol else parameters["cross_parameters"]["Ag_W"]

        if radius <= pair["c_A"]:
            polynomial = sum(pair[f"c{power}"] * radius**power for power in range(5))
            pair_energy += (radius - pair["c_A"]) ** 2 * polynomial

        if radius <= pair["d_A"]:
            delta = radius - pair["d_A"]
            pair_density = pair["A"] ** 2 * (delta**2 + pair["B"] ** 2 * delta**4)
            density[first] += pair_density
            density[second] += pair_density

    return (pair_energy - np.sqrt(density).sum()) / len(atoms)


def run_alloy_crystal(lammps: str, name: str, parameters: dict, work: Path) -> dict:
    directory = work / name
    directory.mkdir()
    a = parameters["lattice_A"]
    atoms = Atoms(
        symbols=parameters["symbols"],
        scaled_positions=parameters["scaled_positions"],
        cell=[a, a, a],
        pbc=True,
    ).repeat((4, 4, 4))
    formula_energy = direct_efs_energy(atoms)
    write(directory / "atoms.data", atoms, format="lammps-data", atom_style="atomic", masses=True, specorder=["Ag", "W"])
    (directory / "in.lmp").write_text(
        f"""units metal
atom_style atomic
boundary p p p
read_data atoms.data
pair_style eam/fs
pair_coeff * * {POTENTIAL} Ag W
neighbor 2.0 bin
thermo 10
thermo_style custom step atoms pe press lx ly lz
run 0
variable initial_epa equal pe/count(all)
print "INITIAL compound {name} epa=$(v_initial_epa:%.12g)"
fix relax all box/relax iso 0.0 vmax 0.001
min_style cg
minimize 1e-12 1e-12 1000 10000
unfix relax
variable epa equal pe/count(all)
variable lat equal lx/4.0
print "RESULT compound {name} epa=$(v_epa:%.12g) lattice=$(v_lat:%.12g)"
"""
    )
    match, output = lammps_result(
        [lammps, "-in", "in.lmp", "-log", "none"],
        directory,
        rf"RESULT compound {name} epa=([-+0-9.eE]+) lattice=([-+0-9.eE]+)",
    )
    cohesive = -float(match.group(1))
    lattice = float(match.group(2))
    initial_match = re.search(rf"INITIAL compound {name} epa=([-+0-9.eE]+)", output)
    if initial_match is None:
        raise RuntimeError(f"LAMMPS did not report the initial {name} energy")
    lammps_initial_energy = float(initial_match.group(1))
    return {
        "cohesive_eV_atom": cohesive,
        "cohesive_table7_this_work_eV_atom": parameters["cohesive_this_work_eV_atom"],
        "cohesive_table7_vasp_eV_atom": parameters["cohesive_vasp_eV_atom"],
        "cohesive_absolute_error_vs_this_work_eV_atom": abs(cohesive - parameters["cohesive_this_work_eV_atom"]),
        "cohesive_absolute_error_vs_vasp_eV_atom": abs(cohesive - parameters["cohesive_vasp_eV_atom"]),
        "direct_efs_initial_energy_eV_atom": formula_energy,
        "lammps_eamfs_initial_energy_eV_atom": lammps_initial_energy,
        "direct_formula_vs_lammps_initial_error_eV_atom": abs(formula_energy - lammps_initial_energy),
        "lattice_A": lattice,
        "lattice_table7_A": parameters["lattice_A"],
        "lattice_absolute_error_A": abs(lattice - parameters["lattice_A"]),
    }


def orient_cell_axes(atoms: Atoms) -> Atoms:
    rotation = (atoms.cell.array / atoms.cell.lengths()[:, None]).T
    atoms.positions = atoms.positions @ rotation
    atoms.set_cell(atoms.cell.array @ rotation)
    return atoms


def build_interface_slab(gap: float) -> tuple[Atoms, float]:
    ag_parameter = REFERENCE["Ag"]["lattice_A"]
    w_parameter = REFERENCE["W"]["lattice_A"]
    ag = orient_cell_axes(make_supercell(
        bulk("Ag", "fcc", a=ag_parameter, cubic=True),
        [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
    )).repeat((2, 2, 3))
    w = orient_cell_axes(make_supercell(
        bulk("W", "bcc", a=w_parameter, cubic=True),
        [[0, 1, 1], [0, -1, 1], [1, 0, 0]],
    )).repeat((2, 2, 3))
    lateral_primitive = (ag_parameter + w_parameter * np.sqrt(2.0)) / 2.0
    lateral_size = 2.0 * lateral_primitive
    ag_positions = ag.positions.copy()
    w_positions = w.positions.copy()
    ag_positions[:, :2] *= lateral_size / ag.cell.lengths()[0]
    w_positions[:, :2] *= lateral_size / w.cell.lengths()[0]
    w_positions[:, :2] += lateral_primitive / 2.0
    w_positions[:, :2] %= lateral_size

    ag_top = float(ag_positions[:, 2].max())
    w_top = float(w_positions[:, 2].max())
    w_positions[:, 2] += ag_top + gap

    box_z = ag_top + gap + w_top + gap
    positions = np.vstack((ag_positions, w_positions))
    cell = Atoms(
        symbols=["Ag"] * len(ag) + ["W"] * len(w),
        positions=positions,
        cell=[lateral_size, lateral_size, box_z],
        pbc=True,
    )
    strain = {
        "Ag": lateral_size / ag.cell.lengths()[0] - 1.0,
        "W": lateral_size / w.cell.lengths()[0] - 1.0,
    }
    cell.info["lateral_strain"] = strain
    return cell, lateral_size**2


def run_interface_gap(lammps: str, gap: float, work: Path) -> dict:
    atoms, area = build_interface_slab(gap)
    directory = work / f"interface_{gap:.1f}"
    directory.mkdir()
    write(
        directory / "atoms.data", atoms, format="lammps-data", atom_style="atomic",
        masses=True, specorder=["Ag", "W"],
    )
    element_symbols = np.asarray(atoms.get_chemical_symbols())
    w_positions = atoms.positions[element_symbols == "W"]
    ag_positions = atoms.positions[element_symbols == "Ag"]
    box = atoms.cell.lengths()
    minimum_distance = float(cKDTree(w_positions, boxsize=box).query(ag_positions)[0].min())
    (directory / "in.lmp").write_text(
        f"""units metal
atom_style atomic
boundary p p p
read_data atoms.data
pair_style eam/fs
pair_coeff * * {POTENTIAL} Ag W
neighbor 2.0 bin
thermo 1
thermo_style custom step atoms pe press
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
run 0
print "RESULT gap={gap:.4f} energy=$(pe:%.12g) fmax=$(c_fmax:%.12g)"
"""
    )
    match, _ = lammps_result(
        [lammps, "-in", "in.lmp", "-log", "none"],
        directory,
        rf"RESULT gap={gap:.4f} energy=([-+0-9.eE]+) fmax=([-+0-9.eE]+)",
    )
    return {
        "gap_A": gap,
        "energy_eV": float(match.group(1)),
        "max_force_eV_A": float(match.group(2)),
        "min_Ag_W_distance_A": minimum_distance,
        "area_A2": area,
        "atoms": len(atoms),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lammps", default=os.environ.get("LAMMPS_EXECUTABLE", "/workspace/.venvs/ag-qsgs/bin/lmp"))
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    lammps = shutil.which(args.lammps) or (args.lammps if Path(args.lammps).is_file() else None)
    if not lammps:
        raise FileNotFoundError(f"LAMMPS executable not found: {args.lammps}")
    if not POTENTIAL.is_file():
        raise FileNotFoundError(f"Generate the potential first: {POTENTIAL}")

    report = {
        "title": "Dai et al. (2006) Ag-W EFS LAMMPS eam/fs validation",
        "potential": str(POTENTIAL.relative_to(ROOT)),
        "potential_sha256": hashlib.sha256(POTENTIAL.read_bytes()).hexdigest(),
        "lammps_version": get_lammps_version(lammps),
        "primary_source_parameter_audit": audit_source_parameters(),
        "pure_crystals": {},
        "cross_compound_benchmarks": {},
        "interface_scan": {},
        "limits": [
            "Cohesive energy is reported as -E_total/N because an isolated atom has zero energy under the paper's EFS formula. Pure-phase reproduction confirms this reference matches the paper's cohesive-energy convention.",
            "An independent analytic evaluation of the Dai equations agrees with the LAMMPS eam/fs energies for B2 and L12 at their reported initial lattice constants.",
            "The small coherent Ag(100)/W(100) scan uses the reported in-plane orientation relationship and a rigid BCC-site registry. Its work of separation is compared with the He et al. DFT target, but it is not a relaxed, converged, independent DFT validation.",
            "The reference pack notes a reported W melting point of 4500 K versus 3695 K experimentally. High-temperature W behavior remains unvalidated.",
            "The AgW_B2 cohesive energy differs from Dai Table 7's EFS result by 0.0859 eV/atom, while differing from its VASP value by 0.0081 eV/atom. The Ag3W_L12 energy differs from the EFS result by 0.0191 eV/atom and from VASP by 0.0060 eV/atom. Resolve the Table 7 discrepancy before production use.",
        ],
    }

    with tempfile.TemporaryDirectory(prefix="dai-efs-validation-") as temporary_directory:
        work = Path(temporary_directory)
        for element in ("Ag", "W"):
            result = run_pure_crystal(lammps, element, work)
            result["pass"] = (
                result["cohesive_absolute_error_eV_atom"] < 0.03
                and result["lattice_absolute_error_A"] < 0.02
            )
            report["pure_crystals"][element] = result
            if not result["pass"]:
                raise AssertionError(f"{element} pure-crystal benchmark failed: {result}")

        for compound, parameters in ALLOY_REFERENCE.items():
            result = run_alloy_crystal(lammps, compound, parameters, work)
            result["close_to_vasp_target"] = (
                result["cohesive_absolute_error_vs_vasp_eV_atom"] < 0.02
                and result["lattice_absolute_error_A"] < 0.02
            )
            result["reproduces_paper_efs_value_within_0p02_eV"] = (
                result["cohesive_absolute_error_vs_this_work_eV_atom"] < 0.02
            )
            report["cross_compound_benchmarks"][compound] = result

        scan = [run_interface_gap(lammps, gap, work) for gap in GAPS_A]
        reference_energy = next(item["energy_eV"] for item in scan if item["gap_A"] == max(GAPS_A))
        area = scan[-1]["area_A2"]
        for item in scan:
            item["binding_energy_eV_A2_vs_6A"] = (item["energy_eV"] - reference_energy) / (2 * area)
            item["work_of_separation_eV_A2_vs_6A"] = -item["binding_energy_eV_A2_vs_6A"]
            item["work_of_separation_J_m2_vs_6A"] = item["work_of_separation_eV_A2_vs_6A"] * 16.021766208
        minimum = min(scan, key=lambda item: item["energy_eV"])
        closest = scan[0]
        report["interface_scan"] = {
            "orientation": "Ag(100)/W(100), [010]Ag parallel [011]W, two periodic interfaces",
            "registry": "BCC-site lateral translation, fixed coordinates",
            "reference_DFT": {
                "citation": "He et al., J. Phys. Chem. C 122 (2018) 17928-17935, Table 4",
                "work_of_separation_J_m2": 2.80,
                "interface_energy_J_m2": 0.72,
            },
            "supercell_lateral_A": build_interface_slab(6.0)[0].cell.lengths()[:2].tolist(),
            "lateral_strain": build_interface_slab(6.0)[0].info["lateral_strain"],
            "fixed_atom_count": 0,
            "gap_reference_A": max(GAPS_A),
            "minimum_energy_gap_A": minimum["gap_A"],
            "minimum_binding_energy_eV_A2": minimum["binding_energy_eV_A2_vs_6A"],
            "maximum_work_of_separation_eV_A2": minimum["work_of_separation_eV_A2_vs_6A"],
            "maximum_work_of_separation_J_m2": minimum["work_of_separation_J_m2_vs_6A"],
            "work_of_separation_error_vs_DFT_J_m2": abs(
                minimum["work_of_separation_J_m2_vs_6A"] - 2.80
            ),
            "shortest_gap_energy_above_minimum": closest["energy_eV"] > minimum["energy_eV"],
            "all_results_finite": all(
                np.isfinite(item["energy_eV"])
                and np.isfinite(item["max_force_eV_A"])
                and item["max_force_eV_A"] >= 0.0
                for item in scan
            ),
            "points": scan,
        }
        if (
            not report["interface_scan"]["all_results_finite"]
            or minimum["gap_A"] == max(GAPS_A)
            or minimum["energy_eV"] >= reference_energy
            or closest["energy_eV"] <= minimum["energy_eV"]
        ):
            raise AssertionError(f"Ag/W interface scan failed: {report['interface_scan']}")

    report["status"] = "PURE_PHASES_AND_ANALYTIC_MAPPING_PASS_TABLE7_B2_DISCREPANCY_REMAINS"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
