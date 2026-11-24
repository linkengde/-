#!/usr/bin/env python3
"""Tiny synthetic smoke test for configurable padded z bounds; no LAMMPS required."""
from __future__ import annotations

import numpy as np

from agw_ttm_reference import map_agw_newton_layer


def main() -> None:
    # One W atom near the lower guard layer and one Ag atom above the old z=200
    # cut. Both must be included by the explicit [35,205] Å mapping interval.
    xyz = np.array([
        [50.0, 50.0, 40.5],
        [100.0, 100.0, 200.5],
    ], dtype=float)
    types = np.array([1, 2], dtype=int)
    mapped = map_agw_newton_layer(
        xyz,
        types,
        grid_shape=(40, 40, 34),
        spacing_A=(5.0, 5.0, 5.0),
        source_kwargs={"allow_above_surface": True},
        newton_mask=np.array([True, True]),
        max_transfer_distance_A=20.0,
        grid_offset_xy_A=(0.0, 0.0),
        box_xy_A=200.0,
        map_z_bounds_A=(35.0, 205.0),
    )
    assert mapped["candidate_newton_atoms"] == 2
    assert len(mapped["atom_indices_ag"]) == 1
    assert len(mapped["atom_indices_w"]) == 1
    assert mapped["excluded_outside_ttm_region"]["below_map_z_min_A"] == 0
    assert mapped["excluded_outside_ttm_region"]["above_map_z_max_A"] == 0
    assert mapped["map_z_bounds_A"] == [35.0, 205.0]
    assert mapped["max_total_fill_fraction"] <= 1.0 + 1e-10
    assert np.isclose(
        mapped["source_power_on_state_eV_ps"],
        float(mapped["source_power_ag_eV_ps"].sum()),
        rtol=0.0,
        atol=1e-10,
    )
    print({
        "status": "PASS_SYNTHETIC_PADDED_Z_SMOKE",
        "mapped_Ag": len(mapped["atom_indices_ag"]),
        "mapped_W": len(mapped["atom_indices_w"]),
        "map_z_bounds_A": mapped["map_z_bounds_A"],
        "max_total_fill": mapped["max_total_fill_fraction"],
        "source_power_eV_ps": mapped["source_power_on_state_eV_ps"],
    })


if __name__ == "__main__":
    main()
