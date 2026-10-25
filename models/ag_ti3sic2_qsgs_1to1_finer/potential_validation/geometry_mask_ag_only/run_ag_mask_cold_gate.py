#!/usr/bin/env python3
"""Run a bounded Ag-only cold gate with an external QSGS geometry wall.

Only Ag--Ag EAM acts atomistically. Ti/Si/C are absent. The QSGS mask is
implemented as a conservative signed-distance penalty through LAMMPS fix
external; wall potential energy is included in the energy ledger.
"""

from __future__ import annotations

import json
import gzip
import csv
import os
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent.parent
sys.path.extend(("/workspace/.venvs/mace/lib/python3.12/site-packages", str(MODEL)))
from lammps import lammps


class SignedDistanceWall:
    def __init__(self, lmp, sdf_path: Path, stiffness: float, mobile_ids: np.ndarray):
        archive = np.load(sdf_path)
        self.phi = archive["phi_A"].astype(np.float64)
        self.spacing = archive["spacing_A"].astype(float)
        self.size = archive["size_A"].astype(float)
        self.clearance = float(archive["clearance_A"])
        self.nx, self.ny, self.nz = self.phi.shape
        self.mobile_tag_mask = np.zeros(int(mobile_ids.max()) + 1, dtype=bool)
        self.mobile_tag_mask[mobile_ids] = True
        self.k = float(stiffness)
        self.lmp = lmp
        self.last = {"active": 0, "max_force_eV_A": 0.0, "energy_eV": 0.0, "max_penetration_A": 0.0}
        self.history = []

    def sample_phi_gradient(self, xyz: np.ndarray):
        """Sample trilinear phi and its exact cell derivative.

        The analytical derivative of the same interpolant used for the energy
        keeps the geometric wall conservative (apart from cell-boundary
        derivative jumps), unlike interpolating a separately tabulated grad.
        """
        dx, dy, dz = self.spacing
        qx = np.mod(xyz[:, 0], self.size[0]) / dx - 0.5
        qy = np.mod(xyz[:, 1], self.size[1]) / dy - 0.5
        qz = np.clip(xyz[:, 2] / dz - 0.5, 0.0, self.nz - 1.0)
        bx, by, bz = np.floor(qx).astype(int), np.floor(qy).astype(int), np.floor(qz).astype(int)
        fx, fy, fz = qx - bx, qy - by, qz - bz
        x0, x1 = bx % self.nx, (bx + 1) % self.nx
        y0, y1 = by % self.ny, (by + 1) % self.ny
        z0, z1 = np.clip(bz, 0, self.nz - 1), np.clip(bz + 1, 0, self.nz - 1)

        v000 = self.phi[x0, y0, z0]; v100 = self.phi[x1, y0, z0]
        v010 = self.phi[x0, y1, z0]; v110 = self.phi[x1, y1, z0]
        v001 = self.phi[x0, y0, z1]; v101 = self.phi[x1, y0, z1]
        v011 = self.phi[x0, y1, z1]; v111 = self.phi[x1, y1, z1]
        wx0, wy0, wz0 = 1.0 - fx, 1.0 - fy, 1.0 - fz
        phi = (v000 * wx0 * wy0 * wz0 + v100 * fx * wy0 * wz0
               + v010 * wx0 * fy * wz0 + v110 * fx * fy * wz0
               + v001 * wx0 * wy0 * fz + v101 * fx * wy0 * fz
               + v011 * wx0 * fy * fz + v111 * fx * fy * fz)
        gx = ((v100 - v000) * wy0 * wz0 + (v110 - v010) * fy * wz0
              + (v101 - v001) * wy0 * fz + (v111 - v011) * fy * fz) / dx
        gy = ((v010 - v000) * wx0 * wz0 + (v110 - v100) * fx * wz0
              + (v011 - v001) * wx0 * fz + (v111 - v101) * fx * fz) / dy
        gz = ((v001 - v000) * wx0 * wy0 + (v101 - v100) * fx * wy0
              + (v011 - v010) * wx0 * fy + (v111 - v110) * fx * fy) / dz
        outside_z = (xyz[:, 2] < 0.0) | (xyz[:, 2] >= self.size[2])
        phi[outside_z] = 1.0e6
        grad = np.column_stack((gx, gy, gz))
        grad[outside_z] = 0.0
        return phi, grad

    def callback(self, _caller, step, nlocal, tag, x, fext):
        xyz = np.asarray(x, dtype=np.float64)
        mobile = self.mobile_tag_mask[np.asarray(tag, dtype=np.int64)]
        phi, gradphi = self.sample_phi_gradient(xyz)
        depth = np.maximum(self.clearance - phi, 0.0)
        depth[~mobile] = 0.0
        active = (depth > 0.0) & mobile
        if np.any(active) and self.k > 0.0:
            force = self.k * depth[:, None] * gradphi
            force[~active] = 0.0
            # LAMMPS reuses the fix-external force buffer between callbacks;
            # overwrite it so this conservative force is not accumulated.
            fext[:] = force
            wall_energy = 0.5 * self.k * float(np.dot(depth, depth))
            self.lmp.fix_external_set_energy_global("geom", wall_energy)
            force_norm = np.linalg.norm(force, axis=1)
            self.last = {
                "active": int(active.sum()),
                "max_force_eV_A": float(force_norm.max(initial=0.0)),
                "energy_eV": float(wall_energy),
                "max_penetration_A": float(depth.max(initial=0.0)),
            }
        else:
            fext[:] = 0.0
            self.lmp.fix_external_set_energy_global("geom", 0.0)
            self.last = {"active": 0, "max_force_eV_A": 0.0, "energy_eV": 0.0, "max_penetration_A": 0.0}
        self.history.append({"step": int(step), **self.last})


def compute_force_metrics(lmp):
    x = np.asarray(lmp.numpy.extract_atom("x"))
    f = np.asarray(lmp.numpy.extract_atom("f"))
    nlocal = int(lmp.extract_setting("nlocal"))
    mobile = x[:nlocal, 2] > 20.0
    fm = np.linalg.norm(f[:nlocal][mobile], axis=1)
    return {
        "mobile_atoms": int(mobile.sum()),
        "force_p50_eV_A": float(np.percentile(fm, 50)),
        "force_p95_eV_A": float(np.percentile(fm, 95)),
        "force_p99_eV_A": float(np.percentile(fm, 99)),
        "force_max_eV_A": float(fm.max()),
    }


def read_dump(path: Path):
    lines = path.read_text().splitlines()
    frame = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS"))
    columns = lines[frame].split()[2:]
    data = np.loadtxt(lines[frame + 1 :], dtype=float)
    result = {name: data[:, i] for i, name in enumerate(columns)}
    order = np.argsort(result["id"].astype(np.int64))
    return {name: val[order] for name, val in result.items()}


def cold_structure_metrics(start_path: Path, end_path: Path, sdf_path: Path):
    start = read_dump(start_path)
    end = read_dump(end_path)
    if not np.array_equal(start["id"], end["id"]):
        raise ValueError("Start/end Ag atom IDs differ")
    xyz0 = np.column_stack([start[c] for c in ("x", "y", "z")])
    xyz1 = np.column_stack([end[c] for c in ("x", "y", "z")])
    delta = xyz1 - xyz0
    box = np.load(sdf_path)
    size = box["size_A"].astype(float)
    for axis in (0, 1):
        delta[:, axis] -= np.rint(delta[:, axis] / size[axis]) * size[axis]
    displacement = np.linalg.norm(delta, axis=1)

    phi_grid = box["phi_A"].astype(float)
    spacing = box["spacing_A"].astype(float)
    def mask_phi(xyz):
        q = np.empty((3, len(xyz)), dtype=float)
        q[0] = np.mod(xyz[:, 0], size[0]) / spacing[0] - 0.5
        q[1] = np.mod(xyz[:, 1], size[1]) / spacing[1] - 0.5
        q[2] = np.clip(xyz[:, 2] / spacing[2] - 0.5, 0.0, phi_grid.shape[2] - 1.0)
        phi = map_coordinates(phi_grid, q, order=1, mode="grid-wrap", prefilter=False)
        phi[(xyz[:, 2] < 0.0) | (xyz[:, 2] >= size[2])] = 1.0e6
        return phi

    phi0, phi1 = mask_phi(xyz0), mask_phi(xyz1)
    boxsize = np.array([size[0], size[1], 1000.0])
    conn = {}
    for label, xyz in (("initial", xyz0), ("final", xyz1)):
        wrapped = xyz.copy()
        wrapped[:, 0] %= size[0]
        wrapped[:, 1] %= size[1]
        tree = cKDTree(wrapped, boxsize=boxsize)
        pairs = tree.query_pairs(3.35, output_type="ndarray")
        n = len(xyz)
        graph = coo_matrix((np.ones(2 * len(pairs), dtype=np.uint8),
                            (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])),
                           shape=(n, n)).tocsr()
        ncomp, labels = connected_components(graph, directed=False, return_labels=True)
        counts = np.bincount(labels)
        conn[label] = {
            "cutoff_A": 3.35,
            "components": int(ncomp),
            "largest_component_atoms": int(counts.max(initial=0)),
            "largest_component_fraction": float(counts.max(initial=0) / n),
            "edges": int(len(pairs)),
        }

    return {
        "atom_count": int(len(xyz0)),
        "displacement_A": {
            "rms": float(np.sqrt(np.mean(displacement**2))),
            "p95": float(np.percentile(displacement, 95)),
            "p99": float(np.percentile(displacement, 99)),
            "max": float(displacement.max(initial=0.0)),
        },
        "atoms_below_z20_A_final": int(np.count_nonzero(xyz1[:, 2] <= 20.0)),
        "mask_signed_distance_A": {
            "initial_min": float(phi0.min()),
            "final_min": float(phi1.min()),
            "initial_inside_count": int(np.count_nonzero(phi0 < 0.0)),
            "final_inside_count": int(np.count_nonzero(phi1 < 0.0)),
        },
        "connectivity": conn,
        "force_magnitudes_eV_A": {
            "initial_p95": float(np.percentile(np.linalg.norm(np.column_stack([start[c] for c in ("fx", "fy", "fz")]), axis=1), 95)),
            "initial_max": float(np.linalg.norm(np.column_stack([start[c] for c in ("fx", "fy", "fz")]), axis=1).max()),
            "final_p95": float(np.percentile(np.linalg.norm(np.column_stack([end[c] for c in ("fx", "fy", "fz")]), axis=1), 95)),
            "final_max": float(np.linalg.norm(np.column_stack([end[c] for c in ("fx", "fy", "fz")]), axis=1).max()),
        },
    }


def thermo(lmp):
    return {
        "step": float(lmp.get_thermo("step")),
        "time": float(lmp.get_thermo("time")),
        "mobile_temp_K": float(lmp.extract_compute("Tmobile", 0, 0)),
        "pe_eV": float(lmp.get_thermo("pe")),
        "ke_eV": float(lmp.get_thermo("ke")),
        "etotal_eV": float(lmp.get_thermo("etotal")),
    }


def summarize_nvt_log(path: Path, steps: int, duration_ps: float):
    lines = path.read_text().splitlines()
    start = max(i for i, line in enumerate(lines) if line.strip() == f"run {steps}")
    rows, reading = [], False
    for line in lines[start:]:
        fields = line.split()
        if fields and fields[0] == "Step" and "c_Tmobile" in line:
            reading = True
            continue
        if reading:
            try:
                rows.append([float(x) for x in fields[:7]])
            except (ValueError, IndexError):
                if rows:
                    break
    data = np.asarray(rows, dtype=float)
    if not len(data):
        raise RuntimeError(f"Could not parse NVT thermodynamic rows from {path}")
    window = data[data[:, 1] >= max(0.0, duration_ps - 0.5)]
    return {
        "sample_count": int(len(data)),
        "final_temperature_K": float(data[-1, 2]),
        "last_0p5ps_temperature_mean_K": float(window[:, 2].mean()),
        "last_0p5ps_temperature_std_K": float(window[:, 2].std()),
        "last_0p5ps_potential_energy_slope_eV_per_ps": float(np.polyfit(window[:, 1], window[:, 3], 1)[0]),
        "last_0p5ps_max_force_peak_eV_A": float(window[:, 6].max()),
    }


def wall_force_finite_difference_check(mask: SignedDistanceWall):
    with gzip.open(HERE / "Ag_only_id_map.csv.gz", "rt", newline="") as f:
        rows = csv.reader(f)
        next(rows)
        xyz = np.asarray([[float(row[2]), float(row[3]), float(row[4])] for row in rows], dtype=float)
    phi, grad = mask.sample_phi_gradient(xyz)
    tags = np.arange(1, len(xyz) + 1)
    selected = np.flatnonzero((phi < mask.clearance) & (xyz[:, 2] > 20.0) & mask.mobile_tag_mask[tags])
    if not len(selected):
        raise RuntimeError("No Ag positions overlap the mask wall for derivative check")
    selected = selected[np.linspace(0, len(selected) - 1, min(128, len(selected)), dtype=int)]
    points = xyz[selected]
    k, eps = 5.0, 1.0e-5
    phi, grad = mask.sample_phi_gradient(points)
    depth = np.maximum(mask.clearance - phi, 0.0)
    analytic_force = k * depth[:, None] * grad
    fd_force = np.zeros_like(analytic_force)
    for axis in range(3):
        plus, minus = points.copy(), points.copy()
        plus[:, axis] += eps
        minus[:, axis] -= eps
        pphi, _ = mask.sample_phi_gradient(plus)
        mphi, _ = mask.sample_phi_gradient(minus)
        ep = 0.5 * k * np.maximum(mask.clearance - pphi, 0.0) ** 2
        em = 0.5 * k * np.maximum(mask.clearance - mphi, 0.0) ** 2
        fd_force[:, axis] = -(ep - em) / (2.0 * eps)
    abs_err = np.abs(fd_force - analytic_force)
    rel_err = abs_err / np.maximum(np.abs(analytic_force), 1.0e-8)
    return {
        "sampled_initial_wall_atoms": int(len(selected)),
        "epsilon_A": eps,
        "max_abs_force_error_eV_A": float(abs_err.max(initial=0.0)),
        "max_relative_force_error": float(rel_err.max(initial=0.0)),
        "rms_force_error_eV_A": float(np.sqrt(np.mean(abs_err**2))),
    }


def load_ag_id_map():
    with gzip.open(HERE / "Ag_only_id_map.csv.gz", "rt", newline="") as f:
        rows = csv.reader(f)
        next(rows)
        data = np.asarray([[int(row[0]), float(row[2]), float(row[3]), float(row[4])] for row in rows])
    tags = data[:, 0].astype(np.int64)
    xyz = data[:, 1:].astype(float)
    return tags, xyz


def main():
    stiffnesses = (0.0, 1.0, 5.0, 10.0)
    fire_mode = "--fire" in sys.argv
    # The selected phase-candidate gets one bounded 2 ps cold screen after FIRE.
    nvt_steps = 4000
    dt_ps = 0.0005
    suffix = "fire_gate_phase_candidate" if fire_mode else "cold_gate_2ps"
    log = HERE / f"log.ag_mask_{suffix}"
    dump_start = HERE / ("ag_mask_phase_fire_nvt_start.dump" if fire_mode else "ag_mask_cold_start.dump")
    dump_end = HERE / ("ag_mask_phase_fire_nvt_final.dump" if fire_mode else "ag_mask_cold_final.dump")

    lmp = lammps(cmdargs=["-screen", "none", "-log", str(log), "-sf", "omp"])
    commands = [
        "units metal",
        "atom_style atomic",
        "boundary p p f",
        "atom_modify sort 0 0",
        f"read_data {HERE}/Ag_only_qsgs.data.gz",
        "pair_style eam",
        f"pair_coeff * * {MODEL}/potentials/Ag_u3.eam",
        "neighbor 2.0 bin",
        "neigh_modify delay 0 every 1 check yes",
        "region r_bottom block INF INF INF INF 0.0 20.0 units box",
        "group bottom region r_bottom",
        "group mobile subtract all bottom",
        "fix freeze_bottom bottom setforce 0.0 0.0 0.0",
        f"fix geom mobile external pf/callback 1 1",
        "fix_modify geom energy yes",
        "timestep 0.0005",
        "compute Tmobile mobile temp",
        "variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)",
        "compute fmax mobile reduce max v_fmag",
        "thermo 20",
        "thermo_style custom step time c_Tmobile pe ke etotal c_fmax",
        "thermo_modify flush yes",
    ]
    for command in commands:
        lmp.command(command)

    tags, parent_xyz = load_ag_id_map()
    mobile_ids = tags[parent_xyz[:, 2] > 20.0]
    mask = SignedDistanceWall(lmp, HERE / "qsgs_signed_distance_A.npz", stiffnesses[0], mobile_ids)
    force_consistency = wall_force_finite_difference_check(mask)
    lmp.set_fix_external_callback("geom", mask.callback, caller=mask)

    static_sensitivity = []
    for k in stiffnesses:
        mask.k = k
        lmp.command("run 0 pre yes post no")
        static_sensitivity.append({
            "stiffness_eV_A2": k,
            "thermo": thermo(lmp),
            "mobile_force": compute_force_metrics(lmp),
            "wall_last_callback": dict(mask.last),
        })

    mask.k = 5.0
    minimization = None
    min_start_dump = HERE / ("ag_mask_phase_fire_min_start.dump" if fire_mode else "ag_mask_fire_min_start.dump")
    min_end_dump = HERE / ("ag_mask_phase_fire_min_end.dump" if fire_mode else "ag_mask_fire_min_end.dump")
    if fire_mode:
        lmp.command(f"write_dump all custom {min_start_dump} id type x y z vx vy vz fx fy fz modify sort id")
        lmp.command("min_style fire")
        lmp.command("min_modify dmax 0.05")
        lmp.command("minimize 0.0 5.0e-2 1000 5000")
        minimization = {
            "thermo": thermo(lmp),
            "mobile_force": compute_force_metrics(lmp),
            "max_iterations": 1000,
            "max_force_evaluations": 5000,
            "energy_tolerance_eV": 0.0,
            "force_tolerance_eV_A": 5.0e-2,
        }
        lmp.command(f"write_dump all custom {min_end_dump} id type x y z vx vy vz fx fy fz modify sort id")
        lmp.command("reset_timestep 0")
    lmp.command("velocity mobile create 300.0 190427 mom yes rot no dist gaussian")
    lmp.command("run 0 pre yes post no")
    initial_thermo = thermo(lmp)
    initial_force = compute_force_metrics(lmp)
    lmp.command(f"write_dump all custom {dump_start} id type x y z vx vy vz fx fy fz modify sort id")
    lmp.command("fix nvt mobile nvt temp 300.0 300.0 0.1")
    lmp.command(f"run {nvt_steps}")
    lmp.command("unfix nvt")
    lmp.command("run 0 pre no post no")
    final_thermo = thermo(lmp)
    final_force = compute_force_metrics(lmp)
    lmp.command(f"write_dump all custom {dump_end} id type x y z vx vy vz fx fy fz modify sort id")
    lmp.close()

    structure = cold_structure_metrics(dump_start, dump_end, HERE / "qsgs_signed_distance_A.npz")
    nvt_summary = summarize_nvt_log(log, nvt_steps, nvt_steps * dt_ps)
    if fire_mode:
        minimization["structure_check"] = cold_structure_metrics(min_start_dump, min_end_dump, HERE / "qsgs_signed_distance_A.npz")
    callback_active = [x["active"] for x in mask.history if x["active"] > 0]
    callback_wall_forces = [x["max_force_eV_A"] for x in mask.history]
    callback_energies = [x["energy_eV"] for x in mask.history]
    result = {
        "status": "AG_ONLY_EAM_QSGS_GEOMETRIC_MASK_BOUNDED_FIRE_COLD_GATE" if fire_mode else "AG_ONLY_EAM_QSGS_GEOMETRIC_MASK_COLD_NVT_GATE",
        "atom_model": "Ag atoms only with Ag_u3.eam; Ti/Si/C excluded from LAMMPS",
        "geometry_constraint": "QSGS ceramic occupancy signed-distance penalty via fix external; no Ag-Ti/Si/C pair terms",
        "fcc_phase_candidate": json.loads((HERE / "fcc_lattice_phase_search.json").read_text())["selected"]
        if (HERE / "fcc_lattice_phase_search.json").exists() else None,
        "dt_fs": dt_ps * 1000.0,
        "nvt_duration_ps": nvt_steps * dt_ps,
        "wall_clearance_A": mask.clearance,
        "wall_stiffness_eV_A2": mask.k,
        "wall_force_finite_difference_check": force_consistency,
        "static_wall_stiffness_sensitivity": static_sensitivity,
        "bounded_fire_minimization": minimization,
        "initial": {"thermo": initial_thermo, "mobile_force": initial_force},
        "final": {"thermo": final_thermo, "mobile_force": final_force},
        "nvt_window_summary": nvt_summary,
        "structure_check": structure,
        "wall_callback_summary": {
            "calls": len(mask.history),
            "calls_with_active_wall": len(callback_active),
            "max_active_atoms_in_call": max(callback_active, default=0),
            "max_wall_force_eV_A": max(callback_wall_forces, default=0.0),
            "max_wall_energy_eV": max(callback_energies, default=0.0),
        },
        "no_heat_or_ttm": True,
        "interpretation_limits": [
            "This finite cold test assesses the Ag-only effective geometric constraint, not a physical Ag-Ti3SiC2 interface.",
            "NVT temperature is controlled; potential-energy drift, force tail, mask penetration, and Ag connectivity must also be assessed.",
            "The NVT run is a finite screening test, not proof of a converged 300 K baseline.",
            "No melting or solidification claim is made by this cold gate.",
        ],
        "outputs": {
            "log": log.name,
            "start_dump": dump_start.name,
            "end_dump": dump_end.name,
            "bounded_fire_min_start_dump": min_start_dump.name if fire_mode else None,
            "bounded_fire_min_end_dump": min_end_dump.name if fire_mode else None,
        },
    }
    metrics_path = HERE / ("ag_mask_phase_fire_gate_metrics.json" if fire_mode else "ag_mask_cold_gate_metrics.json")
    history_path = HERE / ("geometry_wall_callback_history_phase_fire.json" if fire_mode else "geometry_wall_callback_history.json")
    metrics_path.write_text(json.dumps(result, indent=2) + "\n")
    history_path.write_text(json.dumps(mask.history, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
