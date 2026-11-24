#!/usr/bin/env python3
"""Gated LAMMPS bridge for a short Ag-W TTM-MD pilot.

The existing Ag-only LAMMPS/TTM pilot supplied the integration pattern. This
version adds two phase electron grids and refuses to run until a manifest
matches the exact data-file hash. A static PASS is required by default;
EXPLORATORY_ONLY requires an explicit command-line override and is reported as
a software-link check. It is not a production melt-pool solver: the phase-
fraction/CIC map stays fixed during the run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
for parent in HERE.parents:
    if (parent / "models/ag_ti3sic2_qsgs_1to1_finer/ag_ttm_reference.py").is_file():
        PACK = parent
        break
else:  # pragma: no cover - protects use outside this repository
    raise RuntimeError("Ag-W research pack root was not found")

sys.path.insert(0, str(HERE))
from agw_ttm_reference import (  # noqa: E402
    AgWMappingError,
    AgWTTMGrid,
    M_VV2E_METAL,
    map_agw_newton_layer,
    lattice_temperature_by_cell,
    rescale_lattice_kinetic,
)

POTENTIAL = next(
    parent / "potential_files/derived/Dai_2006_EFS_AgW.eam.fs"
    for parent in HERE.parents
    if (parent / "potential_files/derived/Dai_2006_EFS_AgW.eam.fs").is_file()
)
MASS_BY_TYPE = {1: 183.84, 2: 107.8682}  # data-file convention: W=1, Ag=2


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_static_gate(data_path: Path, audit_path: Path, exploratory: bool = False) -> dict:
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    observed_hash = sha256_file(data_path)
    accepted_hash = audit.get("data_sha256")
    if not accepted_hash or accepted_hash.lower() != observed_hash:
        raise SystemExit("Refusing TTM pilot: manifest data_sha256 does not match the supplied data file")
    gate = audit.get("md_start_gate")
    if gate == "PASS":
        pass
    elif exploratory and gate == "EXPLORATORY_ONLY":
        if audit.get("static_audit_status") != "STATIC_EFS_FAIL":
            raise SystemExit("Exploratory override requires the manifest to disclose STATIC_EFS_FAIL")
    else:
        raise SystemExit("Refusing MD: static audit must set md_start_gate to PASS (or use an explicit exploratory manifest and --exploratory)")
    if gate == "PASS":
        force_pass = audit.get("static_force_screen_pass")
        if force_pass is None:
            force_pass = audit.get("forces_eV_A", {}).get("pass")
        if force_pass is not True:
            raise SystemExit("Refusing MD: audit does not record a passed force screen")
        if int(audit.get("atoms_outside_material_z_0_200_A", 0)) != 0:
            raise SystemExit("Refusing MD: audit reports atoms outside the 0-200 A material region")
        ratio = audit.get("mass_ratio_Ag_over_W")
        if ratio is not None and (not math.isfinite(float(ratio)) or abs(float(ratio) - 1.0) > 1.0e-3):
            raise SystemExit("Refusing MD: audited Ag/W mass ratio is outside 1:1 ±0.1%")
    return {"audit": audit, "data_sha256": observed_hash, "gate": gate}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="LAMMPS data file accepted by the static audit")
    parser.add_argument("--audit", type=Path, required=True, help="Matching static audit (PASS) or exploratory-state manifest")
    parser.add_argument("--exploratory", action="store_true", help="Allow an explicit EXPLORATORY_ONLY manifest; results are software-link checks, not accepted MD predictions")
    parser.add_argument("--restart-dump", type=Path, help="Optional LAMMPS custom dump to continue positions and velocities from")
    parser.add_argument("--restart-step", type=int, default=4200, help="Timestep recorded in --restart-dump")
    parser.add_argument("--output", type=Path, required=True, help="Output JSON path")
    parser.add_argument("--pulse-fwhm-ps", type=float, required=True, help="Explicit Gaussian pulse FWHM")
    parser.add_argument("--pulse-center-ps", type=float, required=True, help="Explicit Gaussian pulse center")
    parser.add_argument("--cooling-ps", type=float, required=True, help="Post-pulse cooling interval")
    parser.add_argument("--equilibration-ps", type=float, default=0.0, help="Optional shared 300 K equilibration before mapping; default preserves the supplied restart")
    parser.add_argument("--mapping-only", action="store_true", help="Run LAMMPS setup and the restart-to-grid mapping gate, then exit without a TTM step or laser source")
    parser.add_argument("--grid-A", type=float, default=10.0, help="CIC mesh spacing; box lengths must divide evenly")
    parser.add_argument("--map-z-min-A", type=float, default=40.0, help="Lower z boundary of the TTM mapping mesh")
    parser.add_argument("--map-z-max-A", type=float, default=200.0, help="Upper z boundary of the TTM mapping mesh")
    parser.add_argument("--allow-above-surface", action="store_true", help="Clamp source depth to zero for atoms above the nominal flat source plane; exploratory rough-surface approximation")
    parser.add_argument("--grid-offset-x-A", type=float, default=0.0, help="Periodic x phase of the CIC mesh, in A")
    parser.add_argument("--grid-offset-y-A", type=float, default=0.0, help="Periodic y phase of the CIC mesh, in A")
    parser.add_argument("--box-xy-A", type=float, default=200.0, help="Square periodic x/y box length used by the data/restart and TTM mesh")
    parser.add_argument("--max-map-transfer-A", type=float, default=20.0, help="Maximum joint Ag/W cell-volume transfer distance")
    parser.add_argument("--dt-ps", type=float, default=0.0005, help="LAMMPS metal-unit timestep")
    parser.add_argument("--source-scale", type=float, default=1.0, help="Scale on inherited I0=5 eV/(ps A^2)")
    parser.add_argument("--interface-h-W-m2-K", type=float, default=1.0e9, help="Provisional Ag/W electronic interface conductance")
    parser.add_argument("--seed", type=int, default=20261124)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data_path = args.data.resolve()
    audit_path = args.audit.resolve()
    output_path = args.output.resolve()
    if not data_path.is_file() or not audit_path.is_file():
        raise SystemExit("The requested data file and static audit must both exist")
    if output_path in {data_path, audit_path}:
        raise SystemExit("Output path must not overwrite the supplied structure or its static audit")
    gate = require_static_gate(data_path, audit_path, exploratory=args.exploratory)
    restart_path = args.restart_dump.resolve() if args.restart_dump else None
    if restart_path is not None and not restart_path.is_file():
        raise SystemExit(f"Restart dump not found: {restart_path}")
    if args.exploratory:
        expected_restart = gate["audit"].get("restart_dump_sha256")
        expected_step = gate["audit"].get("restart_dump_timestep")
        if restart_path is None or not expected_restart:
            raise SystemExit("Exploratory TTM mode requires a restart dump whose hash is recorded in the manifest")
        if sha256_file(restart_path).lower() != str(expected_restart).lower():
            raise SystemExit("Exploratory restart dump hash does not match the manifest")
        if expected_step is not None and int(expected_step) != args.restart_step:
            raise SystemExit("Exploratory restart timestep does not match the manifest")
    if args.restart_step < 0:
        raise SystemExit("Restart timestep must be non-negative")
    if not POTENTIAL.is_file():
        raise SystemExit(f"The selected Dai EFS potential is missing: {POTENTIAL}")
    for name, value in (("pulse FWHM", args.pulse_fwhm_ps), ("cooling time", args.cooling_ps),
                        ("grid spacing", args.grid_A), ("box x/y length", args.box_xy_A),
                        ("maximum map transfer distance", args.max_map_transfer_A),
                        ("mapping z minimum", args.map_z_min_A), ("mapping z maximum", args.map_z_max_A),
                        ("timestep", args.dt_ps)):
        if not math.isfinite(value) or value <= 0.0:
            raise SystemExit(f"{name} must be finite and positive")
    if not math.isfinite(args.equilibration_ps) or args.equilibration_ps < 0.0:
        raise SystemExit("Equilibration time must be finite and non-negative")
    if args.mapping_only and args.equilibration_ps != 0.0:
        raise SystemExit("--mapping-only requires zero equilibration so it performs no MD integration")
    if not math.isfinite(args.pulse_center_ps) or args.pulse_center_ps < 0.0:
        raise SystemExit("Pulse center must be finite and non-negative")
    if not math.isfinite(args.source_scale) or args.source_scale < 0.0:
        raise SystemExit("Source scale must be finite and non-negative")
    if not math.isfinite(args.interface_h_W_m2_K) or args.interface_h_W_m2_K < 0.0:
        raise SystemExit("Interface conductance must be finite and non-negative")
    if not math.isfinite(args.grid_offset_x_A) or not math.isfinite(args.grid_offset_y_A):
        raise SystemExit("Grid offsets must be finite")
    grid_offset_xy = (args.grid_offset_x_A % args.box_xy_A, args.grid_offset_y_A % args.box_xy_A)
    nx = int(round(args.box_xy_A / args.grid_A))
    if nx <= 0 or abs(args.box_xy_A / nx - args.grid_A) > args.grid_A * 0.1:
        raise SystemExit("The x/y box must be within 10% of an integer number of nominal grid cells")
    map_depth_A = args.map_z_max_A - args.map_z_min_A
    if map_depth_A <= 0.0 or abs(map_depth_A / args.grid_A - round(map_depth_A / args.grid_A)) > 1.0e-9:
        raise SystemExit("The mapping z interval must be positive and divide evenly by --grid-A")

    # Do not import the LAMMPS bindings until the no-run geometry gate passes.
    try:
        from mpi4py import MPI
    except ImportError:  # a serial LAMMPS Python build does not need mpi4py
        MPI = None
    from lammps import lammps

    if MPI is not None and MPI.COMM_WORLD.Get_size() != 1:
        raise SystemExit("This reference pilot currently requires one MPI rank")
    equil_steps = int(round(args.equilibration_ps / args.dt_ps))
    end_time_ps = args.pulse_center_ps + 4.0 * args.pulse_fwhm_ps + args.cooling_ps
    total_steps = int(math.ceil(end_time_ps / args.dt_ps))
    fwhm2 = args.pulse_fwhm_ps**2
    shape = (nx, nx, int(round(map_depth_A / args.grid_A)))
    spacing = (args.box_xy_A / nx, args.box_xy_A / nx, args.grid_A)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_path = output_path.with_suffix(".lammps.log")
    lmp = lammps(cmdargs=["-screen", "none", "-log", str(log_path)])
    try:
        setup_commands = [
            "units metal",
            "atom_style atomic",
            "boundary p p f",
            "atom_modify sort 0 0",
            f"read_data {data_path}",
            "mass 1 183.84",
            "mass 2 107.8682",
            "pair_style eam/fs",
            f"pair_coeff * * {POTENTIAL} W Ag",
            "neighbor 2.0 bin",
            "neigh_modify delay 0 every 1 check yes",
            f"timestep {args.dt_ps:.12g}",
            # Define static groups from the base data before optional restart
            # coordinates are loaded, preserving original fixed IDs.
            "region fixed_region block INF INF INF INF 0 19.999999 units box",
            "region bath_region block INF INF INF INF 20 39.999999 units box",
            # Integrate atoms that have entered the monitor/vacuum region too;
            # only the 40-200 A material subset is mapped to the TTM mesh.
            "region newton_region block INF INF INF INF 40 300 units box",
            "group fixed region fixed_region",
            "group bath region bath_region",
            "group newton region newton_region",
            "group mobile subtract all fixed",
        ]
        for command in setup_commands:
            lmp.command(command)

        boxlo, boxhi, *_ = lmp.extract_box()
        actual_xy = np.asarray(boxhi[:2], dtype=float) - np.asarray(boxlo[:2], dtype=float)
        if not np.allclose(actual_xy, args.box_xy_A, rtol=0.0, atol=1e-8):
            raise RuntimeError(f"--box-xy-A={args.box_xy_A:g} does not match data box lengths {actual_xy.tolist()}")

        # The LAMMPS groups are assigned from the template data-file positions
        # before read_dump updates coordinates. Preserve that exact static
        # Newton membership for the TTM map and the following integration.
        nlocal = int(lmp.extract_setting("nlocal"))
        natoms = int(round(lmp.get_natoms()))
        if nlocal != natoms:
            raise RuntimeError(f"Expected all {natoms} atoms on one rank, found {nlocal}")
        template_positions_view = lmp.numpy.extract_atom("x")
        template_types_view = lmp.numpy.extract_atom("type")
        template_ids_view = lmp.numpy.extract_atom("id")
        template_positions = np.array(template_positions_view[:nlocal], copy=True)
        template_types = np.array(template_types_view[:nlocal], dtype=int, copy=True)
        template_ids = np.array(template_ids_view[:nlocal], dtype=np.int64, copy=True)
        if np.any(~np.isin(template_types, tuple(MASS_BY_TYPE))):
            raise RuntimeError("The LAMMPS data file must use type 1=W and type 2=Ag")
        newton_template_mask = (template_positions[:, 2] >= 40.0 - 1e-9) & (template_positions[:, 2] <= 300.0 + 1e-9)
        newton_group_ids = template_ids[newton_template_mask]
        template_type_by_id = dict(zip(template_ids.tolist(), template_types.tolist()))

        if restart_path is not None:
            lmp.command(f"read_dump {restart_path} {args.restart_step} x y z vx vy vz box no replace yes")
        continuation_commands = [
            "velocity fixed set 0.0 0.0 0.0",
            "fix hold fixed setforce 0.0 0.0 0.0",
            "thermo 0",
            "thermo_style custom step time atoms temp pe ke etotal",
        ]
        for command in continuation_commands:
            lmp.command(command)

        types_view = lmp.numpy.extract_atom("type")
        ids_view = lmp.numpy.extract_atom("id")
        current_ids = np.array(ids_view[:nlocal], dtype=np.int64, copy=True)
        types = np.array(types_view[:nlocal], dtype=int, copy=True)
        if len(np.unique(current_ids)) != natoms or not np.array_equal(np.sort(current_ids), np.sort(template_ids)):
            raise RuntimeError("read_dump changed the atom-ID set relative to the template data")
        expected_types = np.fromiter((template_type_by_id[int(atom_id)] for atom_id in current_ids), dtype=int, count=nlocal)
        if not np.array_equal(types, expected_types):
            raise RuntimeError("read_dump changed atom types relative to the template data")
        newton_mask = np.isin(current_ids, newton_group_ids)
        lmp.command("run 0 pre no post no")
        force_view = lmp.numpy.extract_atom("f")
        force_vectors = np.array(force_view[:nlocal], dtype=float, copy=True)
        force_norms = np.linalg.norm(force_vectors, axis=1)
        fixed_ids = template_ids[(template_positions[:, 2] >= 0.0) & (template_positions[:, 2] < 20.0)]
        mobile_mask = ~np.isin(current_ids, fixed_ids)
        force_screen = {
            "all_max_eV_A": float(force_norms.max()),
            "mobile_max_eV_A": float(force_norms[mobile_mask].max()),
            "mobile_rms_eV_A": float(np.sqrt(np.mean(force_norms[mobile_mask] ** 2))),
            "mobile_atoms": int(mobile_mask.sum()),
            "not_an_equilibrium_acceptance": True,
        }
        if equil_steps:
            lmp.command("fix eq mobile nvt temp 300.0 300.0 0.1")
            lmp.command(f"run {equil_steps}")
            lmp.command("unfix eq")
        lmp.command("reset_timestep 0")
        lmp.command("fix bath_integrate bath nve")
        lmp.command("fix newton_integrate newton nve")
        lmp.command(f"fix bath_heat bath langevin 300.0 300.0 0.1 {args.seed} tally yes")
        # Keep LAMMPS thermo output off during the Python-coupled loop. The
        # fix tally is read directly after each sampled step; asking LAMMPS's
        # Thermo class to print this scalar in run-0 setup segfaulted in the
        # bundled MPICH build before any TTM step.
        lmp.command("thermo_style custom step time atoms temp pe ke etotal")
        lmp.command("thermo 0")
        # LAMMPS may reallocate atom arrays during the NVT equilibration.
        # Refresh Python/NumPy views before reading or writing atom state.
        positions_view = lmp.numpy.extract_atom("x")
        types_view = lmp.numpy.extract_atom("type")
        ids_view = lmp.numpy.extract_atom("id")
        positions_after_eq = np.array(positions_view[:nlocal], copy=True)
        current_ids = np.array(ids_view[:nlocal], dtype=np.int64, copy=True)
        types = np.array(types_view[:nlocal], dtype=int, copy=True)
        newton_mask = np.isin(current_ids, newton_group_ids)
        try:
            mapped = map_agw_newton_layer(
                positions_after_eq, types, shape, spacing,
                source_kwargs={"allow_above_surface": bool(args.allow_above_surface)}, newton_mask=newton_mask,
                max_transfer_distance_A=args.max_map_transfer_A,
                grid_offset_xy_A=grid_offset_xy,
                box_xy_A=args.box_xy_A,
                map_z_bounds_A=(args.map_z_min_A, args.map_z_max_A),
            )
        except AgWMappingError as exc:
            if not args.mapping_only:
                raise
            result = {
                "status": "AGW_MAPPING_PREFLIGHT_FAIL",
                "purpose": "LAMMPS setup and restart-to-grid mapping gate only; no TTM step or laser source was applied.",
                "static_gate": {"md_start_gate": gate["gate"], "data_sha256": gate["data_sha256"], "audit": str(audit_path), "exploratory_override": bool(args.exploratory)},
                "restart": {"dump": str(restart_path) if restart_path else None, "sha256": sha256_file(restart_path) if restart_path else None, "timestep": args.restart_step},
                "mapping_parameters": {"box_xy_A": args.box_xy_A, "map_z_bounds_A": [args.map_z_min_A, args.map_z_max_A], "grid_shape": list(shape), "spacing_A": list(spacing), "grid_offset_xy_A": list(grid_offset_xy), "max_transfer_distance_A": args.max_map_transfer_A, "allow_above_surface": bool(args.allow_above_surface)},
                "static_Newton_group_atoms": int(newton_mask.sum()),
                "LAMMPS_run0_force_screen": force_screen,
                "current_Newton_group_atoms_outside_map_z_bounds": int(np.count_nonzero(newton_mask & ((positions_after_eq[:, 2] < args.map_z_min_A - 1e-9) | (positions_after_eq[:, 2] > args.map_z_max_A + 1e-9)))),
                "reason": str(exc),
                "ttm_steps_executed": 0,
                "laser_source_applied": False,
            }
            output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
            print(json.dumps(result, indent=2))
            return
        if args.mapping_only:
            result = {
                "status": "AGW_MAPPING_PREFLIGHT_PASS",
                "purpose": "LAMMPS setup and restart-to-grid mapping gate only; no TTM step or laser source was applied.",
                "static_gate": {"md_start_gate": gate["gate"], "data_sha256": gate["data_sha256"], "audit": str(audit_path), "exploratory_override": bool(args.exploratory)},
                "restart": {"dump": str(restart_path) if restart_path else None, "sha256": sha256_file(restart_path) if restart_path else None, "timestep": args.restart_step},
                "mapping_parameters": {"box_xy_A": args.box_xy_A, "map_z_bounds_A": [args.map_z_min_A, args.map_z_max_A], "grid_shape": list(shape), "spacing_A": list(spacing), "grid_offset_xy_A": list(grid_offset_xy), "max_transfer_distance_A": args.max_map_transfer_A, "allow_above_surface": bool(args.allow_above_surface)},
                "LAMMPS_run0_force_screen": force_screen,
                "mapped_geometry": {"candidate_newton_atoms": mapped["candidate_newton_atoms"], "Ag_atoms": int(len(mapped["atom_indices_ag"])), "W_atoms": int(len(mapped["atom_indices_w"])), "map_z_bounds_A": mapped["map_z_bounds_A"], "excluded_outside_ttm_region": mapped["excluded_outside_ttm_region"], "Ag_source_quadrature_S1_eV_ps": mapped["source_power_on_state_eV_ps"], "max_total_phase_fill": mapped["max_total_fill_fraction"], "joint_mapping_diagnostics": mapped["mapping_diagnostics"]},
                "ttm_steps_executed": 0,
                "laser_source_applied": False,
            }
            output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
            print(json.dumps(result, indent=2))
            return
        idx_ag = mapped["atom_indices_ag"]
        idx_w = mapped["atom_indices_w"]
        masses_all = np.where(types == 1, MASS_BY_TYPE[1], MASS_BY_TYPE[2])
        masses_ag, masses_w = masses_all[idx_ag], masses_all[idx_w]
        electrons = AgWTTMGrid(
            mapped["fraction_ag"], mapped["fraction_w"], spacing,
            source_power_ag_eV_ps=mapped["source_power_ag_eV_ps"],
            interface_conductance_W_m2_K=args.interface_h_W_m2_K,
            initial_temperature_K=300.0,
        )
        atom_map_ag, atom_map_w = mapped["atom_map_ag"], mapped["atom_map_w"]
        velocity_view = lmp.numpy.extract_atom("v")
        velocity = np.array(velocity_view[:nlocal], copy=True)
        atom_energy_initial = float(lmp.get_thermo("pe")) + 0.5 * M_VV2E_METAL * float(
            np.sum(masses_all[:, None] * velocity**2)
        )
        electron_energy_initial = float(
            electrons.electron_energy_ag_eV.sum() + electrons.electron_energy_w_eV.sum()
        )
        bath_energy_start = float(lmp.get_thermo("f_bath_heat"))
        source_total = lattice_ag_total = lattice_w_total = 0.0
        grid_residual_max = 0.0
        combined_residual_max = 0.0
        history: list[dict] = []
        report_every = max(1, total_steps // 100)

        for step in range(1, total_steps + 1):
            lmp.command("run 1 pre no post no")
            velocity_view = lmp.numpy.extract_atom("v")
            velocity = np.array(velocity_view[:nlocal], copy=True)
            tl_ag = lattice_temperature_by_cell(velocity[idx_ag], masses_ag, atom_map_ag, 1)
            tl_w = lattice_temperature_by_cell(velocity[idx_w], masses_w, atom_map_w, 1)
            time_mid = (step - 0.5) * args.dt_ps
            envelope = math.exp(-4.0 * math.log(2.0) * ((time_mid - args.pulse_center_ps) ** 2) / fwhm2)
            electrons.set_source_envelope(args.source_scale * envelope)
            account = electrons.advance(args.dt_ps, tl_ag, tl_w, heat_on=True)
            velocity_ag, applied_ag = rescale_lattice_kinetic(
                velocity[idx_ag], masses_ag, atom_map_ag, account.lattice_energy_ag_eV
            )
            velocity_w, applied_w = rescale_lattice_kinetic(
                velocity[idx_w], masses_w, atom_map_w, account.lattice_energy_w_eV
            )
            velocity[idx_ag] = velocity_ag
            velocity[idx_w] = velocity_w
            velocity_view[:nlocal] = velocity
            source_total += account.source_to_electron_eV
            lattice_ag_total += float(applied_ag.sum())
            lattice_w_total += float(applied_w.sum())
            grid_residual_max = max(grid_residual_max, abs(account.conservation_residual_eV))

            if step == 1 or step % report_every == 0 or step == total_steps:
                pe = float(lmp.get_thermo("pe"))
                kinetic = 0.5 * M_VV2E_METAL * float(np.sum(masses_all[:, None] * velocity**2))
                atom_energy = pe + kinetic
                electron_energy = float(
                    electrons.electron_energy_ag_eV.sum() + electrons.electron_energy_w_eV.sum()
                )
                bath_delta = float(lmp.get_thermo("f_bath_heat")) - bath_energy_start
                residual = atom_energy + electron_energy - atom_energy_initial - electron_energy_initial - source_total - bath_delta
                combined_residual_max = max(combined_residual_max, abs(residual))
                history.append({
                    "step": step, "time_ps": step * args.dt_ps,
                    "pulse_envelope_midpoint": envelope,
                    "max_Te_Ag_K": float(electrons.temperature_ag_K.max()),
                    "max_Te_W_K": float(electrons.temperature_w_K.max()),
                    "source_energy_eV": source_total,
                    "electron_to_lattice_Ag_eV": lattice_ag_total,
                    "electron_to_lattice_W_eV": lattice_w_total,
                    "bath_energy_added_eV": bath_delta,
                    "combined_energy_residual_eV": residual,
                })

        final = history[-1]
        energy_scale_eV = abs(source_total) + abs(final["bath_energy_added_eV"])
        energy_tolerance_eV = max(1.0e-3, 1.0e-3 * energy_scale_eV)
        combined_ledger_pass = combined_residual_max <= energy_tolerance_eV
        grid_ledger_pass = grid_residual_max <= max(1.0e-9, 1.0e-10 * max(abs(source_total), 1.0))
        energy_status = "PASS" if combined_ledger_pass and grid_ledger_pass else "FAIL"
        status_prefix = "AGW_TTM_SOFTWARE_LINK_" if args.exploratory else "AGW_TTM_MD_PILOT_"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result = {
            "status": status_prefix + energy_status,
            "purpose": "Short exploratory LAMMPS/TTM coupling run; not a validated melt/quench prediction.",
            "static_gate": {"md_start_gate": gate["gate"], "data_sha256": gate["data_sha256"], "audit": str(audit_path), "exploratory_override": bool(args.exploratory)},
            "restart": {"dump": str(restart_path) if restart_path else None, "sha256": sha256_file(restart_path) if restart_path else None},
            "atom_model": "Dai 2006 EFS eam/fs; type 1=W, type 2=Ag; p p f; fixed 0-20 A; 300 K bath 20-40 A; TTM/material map 40-200 A; Newton integration 40-300 A to keep monitor-region atoms dynamic",
            "electron_model": "Reused Ag finite-volume reference core with separate Ag/W phase grids and provisional equal-opposite interface exchange",
            "source_model": "Ag-only Gaussian in x/y and exponential depth; temporal Gaussian with requested FWHM and center",
            "spatial_source_parameters": {"I0_eV_ps_A2": 5.0 * args.source_scale, "R_A": 30.0, "skin_depth_A": 20.0, "center_xy_A": [100.0, 100.0], "surface_z_A": 200.0},
            "temporal_parameters": {"pulse_fwhm_ps": args.pulse_fwhm_ps, "pulse_center_ps": args.pulse_center_ps, "cooling_after_center_plus_4fwhm_ps": args.cooling_ps},
            "numerics": {"dt_ps": args.dt_ps, "equilibration_steps": equil_steps, "heating_and_cooling_steps": total_steps, "box_xy_A": args.box_xy_A, "grid_shape": list(shape), "spacing_A": list(spacing), "grid_offset_xy_A": list(grid_offset_xy), "interface_h_W_m2_K": args.interface_h_W_m2_K, "mpi4py_available": MPI is not None, "mpi_ranks": int(MPI.COMM_WORLD.Get_size()) if MPI is not None else 1},
            "mapped_geometry": {"Ag_atoms_Newton": int(len(idx_ag)), "W_atoms_Newton": int(len(idx_w)), "candidate_newton_atoms": mapped["candidate_newton_atoms"], "excluded_outside_ttm_region": mapped["excluded_outside_ttm_region"], "Ag_source_on_state_eV_ps": mapped["source_power_on_state_eV_ps"], "max_total_phase_fill": mapped["max_total_fill_fraction"], "joint_mapping_diagnostics": mapped["mapping_diagnostics"]},
            "energy_check": {"status": energy_status, "source_eV": source_total, "bath_energy_added_eV": final["bath_energy_added_eV"], "final_combined_residual_eV": final["combined_energy_residual_eV"], "max_sampled_combined_residual_eV": combined_residual_max, "combined_tolerance_eV": energy_tolerance_eV, "max_grid_step_residual_eV": grid_residual_max, "combined_ledger_pass": combined_ledger_pass, "grid_ledger_pass": grid_ledger_pass},
            "history": history,
            "limitations": [
                "This pilot run used an exploratory-only override because the v5 EFS static audit remains STATIC_EFS_FAIL; an energy-ledger PASS is a software/numerical coupling result, not a static-structure or physical-validity pass.",
                "The Ag/W phase-fraction and atom-to-grid maps are frozen at the post-read initial coordinates; large melt migration or composition redistribution is not represented.",
                "The Ag/W electronic interface conductance is a provisional analogue-based 1 GW m^-2 K^-1 unless explicitly changed; it is not a measured Ag/W value.",
                "The deterministic velocity-rescaling handoff conserves kinetic energy but does not preserve local momentum or add stochastic phonon fluctuations.",
                "The combined ledger uses the Langevin fix energy tally; this prototype has no matched no-TTM NVE drift control run.",
                "The pulse is a user-specified exploratory Gaussian and is not normalized to a measured absorbed fluence.",
                "A completed pilot is not by itself evidence that the EFS potential quantitatively predicts Ag-W melting or solidification.",
            ],
        }
        output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({key: result[key] for key in ("status", "static_gate", "mapped_geometry", "energy_check", "limitations")}, indent=2))
        if energy_status != "PASS":
            raise SystemExit("Ag-W TTM-MD pilot energy ledger failed the recorded numerical tolerance")
    finally:
        lmp.close()


if __name__ == "__main__":
    main()
