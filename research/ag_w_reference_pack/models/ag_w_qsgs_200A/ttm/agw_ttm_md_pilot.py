#!/usr/bin/env python3
"""Gated LAMMPS bridge for a short Ag-W TTM-MD pilot.

The existing Ag-only LAMMPS/TTM pilot supplied the integration pattern. This
version adds two phase electron grids and refuses to run until a static audit
explicitly accepts the exact data-file hash. It is not a production melt-pool
solver: the phase-fraction/CIC map stays fixed during the run.
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
    AgWTTMGrid,
    M_VV2E_METAL,
    map_agw_newton_layer,
    lattice_temperature_by_cell,
    rescale_lattice_kinetic,
)

POTENTIAL = PACK / "potential_files/derived/Dai_2006_EFS_AgW.eam.fs"
MASS_BY_TYPE = {1: 183.84, 2: 107.8682}  # data-file convention: W=1, Ag=2


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_static_gate(data_path: Path, audit_path: Path) -> dict:
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("md_start_gate") != "PASS":
        raise SystemExit("Refusing MD: static audit must set md_start_gate to PASS")
    observed_hash = sha256_file(data_path)
    accepted_hash = audit.get("data_sha256")
    if not accepted_hash or accepted_hash.lower() != observed_hash:
        raise SystemExit("Refusing MD: static audit data_sha256 does not match the supplied data file")
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
    return {"audit": audit, "data_sha256": observed_hash}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="LAMMPS data file accepted by the static audit")
    parser.add_argument("--audit", type=Path, required=True, help="Static audit JSON with md_start_gate=PASS and data_sha256")
    parser.add_argument("--output", type=Path, required=True, help="Output JSON path")
    parser.add_argument("--pulse-fwhm-ps", type=float, required=True, help="Explicit Gaussian pulse FWHM")
    parser.add_argument("--pulse-center-ps", type=float, required=True, help="Explicit Gaussian pulse center")
    parser.add_argument("--cooling-ps", type=float, required=True, help="Post-pulse cooling interval")
    parser.add_argument("--equilibration-ps", type=float, default=1.0, help="300 K mobile-layer equilibration time")
    parser.add_argument("--grid-A", type=float, default=10.0, help="CIC mesh spacing; box lengths must divide evenly")
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
    gate = require_static_gate(data_path, audit_path)
    if not POTENTIAL.is_file():
        raise SystemExit(f"The selected Dai EFS potential is missing: {POTENTIAL}")
    for name, value in (("pulse FWHM", args.pulse_fwhm_ps), ("cooling time", args.cooling_ps),
                        ("equilibration time", args.equilibration_ps), ("grid spacing", args.grid_A),
                        ("timestep", args.dt_ps)):
        if not math.isfinite(value) or value <= 0.0:
            raise SystemExit(f"{name} must be finite and positive")
    if not math.isfinite(args.pulse_center_ps) or args.pulse_center_ps < 0.0:
        raise SystemExit("Pulse center must be finite and non-negative")
    if not math.isfinite(args.source_scale) or args.source_scale < 0.0:
        raise SystemExit("Source scale must be finite and non-negative")
    if not math.isfinite(args.interface_h_W_m2_K) or args.interface_h_W_m2_K < 0.0:
        raise SystemExit("Interface conductance must be finite and non-negative")
    if args.grid_A <= 0.0 or any(abs(length / args.grid_A - round(length / args.grid_A)) > 1.0e-9
                                 for length in (200.0, 200.0, 160.0)):
        raise SystemExit("Grid spacing must divide the 200x200x160 A Newton-layer mesh exactly")

    # Do not import the LAMMPS bindings until the no-run geometry gate passes.
    from mpi4py import MPI
    from lammps import lammps

    if MPI.COMM_WORLD.Get_size() != 1:
        raise SystemExit("This reference pilot currently requires one MPI rank")
    equil_steps = int(round(args.equilibration_ps / args.dt_ps))
    end_time_ps = args.pulse_center_ps + 4.0 * args.pulse_fwhm_ps + args.cooling_ps
    total_steps = int(math.ceil(end_time_ps / args.dt_ps))
    fwhm2 = args.pulse_fwhm_ps**2
    shape = (int(round(200.0 / args.grid_A)), int(round(200.0 / args.grid_A)), int(round(160.0 / args.grid_A)))
    spacing = (args.grid_A, args.grid_A, args.grid_A)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_path = output_path.with_suffix(".lammps.log")
    lmp = lammps(cmdargs=["-screen", "none", "-log", str(log_path)])
    try:
        commands = [
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
            "region fixed_region block INF INF INF INF 0 19.999999 units box",
            "region bath_region block INF INF INF INF 20 39.999999 units box",
            "region newton_region block INF INF INF INF 40 200 units box",
            "group fixed region fixed_region",
            "group bath region bath_region",
            "group newton region newton_region",
            "group mobile subtract all fixed",
            "velocity fixed set 0.0 0.0 0.0",
            "fix hold fixed setforce 0.0 0.0 0.0",
            "thermo 1",
            "thermo_style custom step time atoms temp pe ke etotal",
        ]
        for command in commands:
            lmp.command(command)

        nlocal = int(lmp.extract_setting("nlocal"))
        natoms = int(round(lmp.get_thermo("atoms")))
        if nlocal != natoms:
            raise RuntimeError(f"Expected all {natoms} atoms on one rank, found {nlocal}")
        positions_view = lmp.numpy.extract_atom("x")
        types_view = lmp.numpy.extract_atom("type")
        positions_initial = np.array(positions_view[:nlocal], copy=True)
        types = np.array(types_view[:nlocal], dtype=int, copy=True)
        if np.any(~np.isin(types, tuple(MASS_BY_TYPE))):
            raise RuntimeError("The LAMMPS data file must use type 1=W and type 2=Ag")
        # Keep membership consistent with the statically defined LAMMPS
        # Newton group during the subsequent 300 K equilibration.
        newton_mask = (positions_initial[:, 2] >= 40.0 - 1e-9) & (positions_initial[:, 2] <= 200.0 + 1e-9)
        lmp.command("run 0 pre no post no")
        lmp.command("fix eq mobile nvt temp 300.0 300.0 0.1")
        lmp.command(f"run {equil_steps}")
        lmp.command("unfix eq")
        lmp.command("reset_timestep 0")
        lmp.command("fix bath_integrate bath nve")
        lmp.command("fix newton_integrate newton nve")
        lmp.command(f"fix bath_heat bath langevin 300.0 300.0 0.1 {args.seed} tally yes")
        lmp.command("thermo_style custom step time atoms temp pe ke etotal f_bath_heat")
        lmp.command("run 0 pre no post no")
        positions_after_eq = np.array(positions_view[:nlocal], copy=True)
        mapped = map_agw_newton_layer(
            positions_after_eq, types, shape, spacing,
            source_kwargs={"allow_above_surface": True}, newton_mask=newton_mask,
        )
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
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result = {
            "status": "AGW_TTM_MD_PILOT_" + energy_status,
            "purpose": "Short exploratory LAMMPS/TTM coupling run; not a validated melt/quench prediction.",
            "static_gate": {"md_start_gate": "PASS", "data_sha256": gate["data_sha256"], "audit": str(audit_path)},
            "atom_model": "Dai 2006 EFS eam/fs; type 1=W, type 2=Ag; p p f; fixed 0-20 A; 300 K bath 20-40 A; Newton 40-200 A",
            "electron_model": "Reused Ag finite-volume reference core with separate Ag/W phase grids and provisional equal-opposite interface exchange",
            "source_model": "Ag-only Gaussian in x/y and exponential depth; temporal Gaussian with requested FWHM and center",
            "spatial_source_parameters": {"I0_eV_ps_A2": 5.0 * args.source_scale, "R_A": 30.0, "skin_depth_A": 20.0, "center_xy_A": [100.0, 100.0], "surface_z_A": 200.0},
            "temporal_parameters": {"pulse_fwhm_ps": args.pulse_fwhm_ps, "pulse_center_ps": args.pulse_center_ps, "cooling_after_center_plus_4fwhm_ps": args.cooling_ps},
            "numerics": {"dt_ps": args.dt_ps, "equilibration_steps": equil_steps, "heating_and_cooling_steps": total_steps, "grid_shape": list(shape), "spacing_A": list(spacing), "interface_h_W_m2_K": args.interface_h_W_m2_K},
            "mapped_geometry": {"Ag_atoms_Newton": int(len(idx_ag)), "W_atoms_Newton": int(len(idx_w)), "Ag_source_on_state_eV_ps": mapped["source_power_on_state_eV_ps"], "max_total_phase_fill": mapped["max_total_fill_fraction"]},
            "energy_check": {"status": energy_status, "source_eV": source_total, "bath_energy_added_eV": final["bath_energy_added_eV"], "final_combined_residual_eV": final["combined_energy_residual_eV"], "max_sampled_combined_residual_eV": combined_residual_max, "combined_tolerance_eV": energy_tolerance_eV, "max_grid_step_residual_eV": grid_residual_max, "combined_ledger_pass": combined_ledger_pass, "grid_ledger_pass": grid_ledger_pass},
            "history": history,
            "limitations": [
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
