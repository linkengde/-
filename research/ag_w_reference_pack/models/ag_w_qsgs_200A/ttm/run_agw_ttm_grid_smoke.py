#!/usr/bin/env python3
"""Small synthetic two-phase grid check; not an Ag–W laser/MD run."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from agw_ttm_reference import (
    AG_ATOMIC_VOLUME_A3,
    AgWTTMGrid,
    BaselineMetalParameters,
    gaussian_ag_source_power,
    map_agw_newton_layer,
)


HERE = Path(__file__).resolve().parent
OUT = HERE / "agw_ttm_grid_smoke_metrics.json"


def main() -> None:
    shape = (8, 8, 4)
    spacing = (10.0, 10.0, 10.0)
    f_ag = np.zeros(shape)
    f_w = np.zeros(shape)
    f_ag[:4, :, :] = 1.0
    f_w[4:, :, :] = 1.0
    source = np.zeros(shape)
    source[3, 4, 2] = 0.5  # eV/ps to one Ag cell; deliberately reduced test input

    ag300 = BaselineMetalParameters("Ag")
    w300 = BaselineMetalParameters("W")
    assert np.isclose(ag300.ke(300.0) / 6.241509074e-4, 429.0, rtol=1e-12)
    assert np.isclose(w300.ke(300.0) / 6.241509074e-4, 130.5, rtol=1e-12)

    grid = AgWTTMGrid(
        f_ag,
        f_w,
        spacing,
        source_power_ag_eV_ps=source,
        interface_conductance_W_m2_K=1.0e9,
        initial_temperature_K=300.0,
    )

    synthetic_atoms = np.array([[100.0, 100.0, 199.0], [110.0, 110.0, 195.0]])
    synthetic_types = np.array([2, 1])  # type 2 Ag, type 1 W
    mapped = map_agw_newton_layer(
        synthetic_atoms,
        synthetic_types,
        (20, 20, 16),
        (10.0, 10.0, 10.0),
    )
    expected_on_state_power = float(
        gaussian_ag_source_power(synthetic_atoms[:1])[0] * AG_ATOMIC_VOLUME_A3
    )
    if not np.isclose(mapped["source_power_on_state_eV_ps"], expected_on_state_power, rtol=1e-12):
        raise RuntimeError("Ag-only source mapping failed its conservative quadrature check")
    e0 = float(grid.electron_energy_ag_eV.sum() + grid.electron_energy_w_eV.sum())
    lattice_reservoir_eV = 0.0
    source_total_eV = 0.0
    interface_total_eV = 0.0
    max_step_residual_eV = 0.0
    max_substeps = 0
    tl_ag = np.full(shape, np.nan)
    tl_w = np.full(shape, np.nan)
    tl_ag[f_ag > 0.0] = 300.0
    tl_w[f_w > 0.0] = 300.0
    steps = 100
    dt_ps = 0.0005
    for step in range(steps):
        grid.set_source_envelope(1.0 if step < 20 else 0.0)
        account = grid.advance(
            dt_ps,
            tl_ag,
            tl_w,
            heat_on=True,
        )
        lattice_reservoir_eV += account.electron_to_lattice_eV
        source_total_eV += account.source_to_electron_eV
        interface_total_eV += account.interface_net_eV
        max_step_residual_eV = max(max_step_residual_eV, abs(account.conservation_residual_eV))
        max_substeps = max(max_substeps, account.substeps)

    e1 = float(grid.electron_energy_ag_eV.sum() + grid.electron_energy_w_eV.sum())
    residual = e1 + lattice_reservoir_eV - e0 - source_total_eV
    result = {
        "status": "PASS" if abs(residual) < 1e-9 and abs(interface_total_eV) < 1e-12 else "FAIL",
        "purpose": "Parameter-load, Ag-only source, two-phase interface exchange, and energy-ledger smoke check; not a physical pulse or MD trajectory.",
        "mesh_shape": list(shape),
        "spacing_A": list(spacing),
        "interface_conductance_W_m2_K": 1.0e9,
        "interface_model": "face-based unlike-phase transfer; equal-and-opposite; no sub-cell interface term",
        "inactive_phase_lattice_temperature_input": "NaN in empty-phase cells; finite values in active cells",
        "test_source": "0.5 eV/ps to one Ag cell for 0.01 ps; then source off for 0.04 ps",
        "dt_ps": dt_ps,
        "steps": steps,
        "source_energy_eV": source_total_eV,
        "electron_to_lattice_energy_eV": lattice_reservoir_eV,
        "initial_electron_energy_eV": e0,
        "final_electron_energy_eV": e1,
        "interface_net_energy_eV": interface_total_eV,
        "synthetic_Ag_only_source_mapping_power_eV_ps": mapped["source_power_on_state_eV_ps"],
        "synthetic_Ag_only_source_expected_power_eV_ps": expected_on_state_power,
        "synthetic_mapping_max_total_phase_fill": mapped["max_total_fill_fraction"],
        "combined_energy_residual_eV": residual,
        "max_step_internal_residual_eV": max_step_residual_eV,
        "max_interface_substeps_per_step": max_substeps,
        "final_max_Te_Ag_K": float(grid.temperature_ag_K.max()),
        "final_max_Te_W_K": float(grid.temperature_w_K.max()),
        "room_temperature_ke_Ag_W_m_K": float(ag300.ke(300.0) / 6.241509074e-4),
        "room_temperature_ke_W_W_m_K": float(w300.ke(300.0) / 6.241509074e-4),
        "reused_core": "models/ag_ti3sic2_qsgs_1to1_finer/ag_ttm_reference.py",
        "limitations": [
            "No atoms or LAMMPS steps are used by this grid-only check.",
            "The test pulse is reduced and numerical only; it does not set the physical pulse envelope or absorbed fluence.",
            "Ag/W interface conductance is a provisional 1 GW m^-2 K^-1 sensitivity value, not an Ag/W measurement.",
            "Full-size QSGS static acceptance and TTM-MD coupling remain outstanding.",
        ],
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
