"""Audit the atom-count and thermal-mask Ag–Ti3SiC2 composition conventions."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import numpy as np


RUN_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUN_DIR.parents[1]
THERMAL_DIR = PROJECT_DIR / "phase_field_qsgs"
TARGET_GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
TARGET_GEOMETRY_REPORT = RUN_DIR / "composition_70vol_geometry.json"
MESH_COMPARISON = RUN_DIR / "mesh2_ag70_comparison.json"


def main() -> None:
    geometry_path = PROJECT_DIR / "qsgs_morphology.npz"
    with np.load(geometry_path, allow_pickle=False) as archive:
        skeleton = np.asarray(archive["skeleton"], dtype=bool)
        domain_size_A = np.asarray(archive["size_A"], dtype=float)
    properties = json.loads((THERMAL_DIR / "thermophysical_candidates.json").read_text())
    atomic = json.loads((PROJECT_DIR / "independent_geometry_validation.json").read_text())
    target_geometry = json.loads(TARGET_GEOMETRY_REPORT.read_text())
    mesh_comparison = json.loads(MESH_COMPARISON.read_text())

    ag_fraction = float(np.mean(~skeleton))
    ceramic_fraction = float(np.mean(skeleton))
    rho_ag = float(properties["ag"]["solid_reference_density_kg_m3"])
    rho_ceramic = float(properties["ti3sic2"]["density_kg_m3"])
    implied_mass_ag = ag_fraction * rho_ag
    implied_mass_ceramic = ceramic_fraction * rho_ceramic
    implied_mass_ratio = implied_mass_ag / implied_mass_ceramic
    ag_volume_fraction_for_equal_bulk_density_mass = rho_ceramic / (rho_ag + rho_ceramic)
    user_ag_volume_fraction = 0.70
    user_ceramic_volume_fraction = 0.30
    user_volume_mass_ratio = (
        user_ag_volume_fraction * rho_ag / (user_ceramic_volume_fraction * rho_ceramic)
    )

    report = {
        "status": "VOLUME_TARGET_CONFIRMED_MASS_CONVERSION_CORRECTED",
        "geometry_path": str(geometry_path.relative_to(PROJECT_DIR)),
        "geometry_sha256": hashlib.sha256(geometry_path.read_bytes()).hexdigest(),
        "mask_shape_xyz": list(skeleton.shape),
        "domain_size_A": domain_size_A.tolist(),
        "thermal_mask_volume_fraction": {
            "Ag": ag_fraction,
            "Ti3SiC2": ceramic_fraction,
        },
        "candidate_bulk_densities_kg_m3": {
            "Ag_solid_reference": rho_ag,
            "Ti3SiC2": rho_ceramic,
        },
        "bulk_density_mass_ratio_implied_by_thermal_mask_Ag_over_Ti3SiC2": implied_mass_ratio,
        "implied_Ag_mass_fraction_by_bulk_density": implied_mass_ag / (implied_mass_ag + implied_mass_ceramic),
        "Ag_volume_fraction_for_equal_mass_at_candidate_bulk_densities": ag_volume_fraction_for_equal_bulk_density_mass,
        "user_stated_target_volume_fraction": {
            "Ag": user_ag_volume_fraction,
            "Ti3SiC2": user_ceramic_volume_fraction,
        },
        "mass_ratio_implied_by_user_stated_70_30_volume_fraction": user_volume_mass_ratio,
        "Ag_mass_fraction_implied_by_user_stated_70_30_volume_fraction": (
            user_ag_volume_fraction * rho_ag
            / (user_ag_volume_fraction * rho_ag + user_ceramic_volume_fraction * rho_ceramic)
        ),
        "atomistic_geometry_mass_ratio_Ag_over_Ti3SiC2": float(atomic["mass_ratio_Ag_over_Ti3SiC2"]),
        "atomistic_geometry_counts": atomic["counts"],
        "interpretation": (
            "The user confirmed Ag:Ti3SiC2 = 70:30 by volume as the experimental target. Under the current "
            "candidate bulk densities this corresponds to about 5.42:1 by mass, not 1:1. The existing "
            "atomistic geometry's approximately equal-mass accounting is a separate construction and does not "
            "set the chosen continuum volume target. The exact density/porosity basis remains to be sourced."
        ),
        "experimental_target_Ag_volume_fraction": user_ag_volume_fraction,
        "user_confirmed_target_basis": "volume_fraction",
        "user_confirmed_Ag_volume_fraction": user_ag_volume_fraction,
        "user_confirmed_Ti3SiC2_volume_fraction": user_ceramic_volume_fraction,
        "claimed_equal_mass_conversion_matches_candidate_density_conversion": False,
        "mass_ratio_from_candidate_density_for_selected_volume_target_Ag_over_Ti3SiC2": user_volume_mass_ratio,
        "selected_target_geometry": target_geometry,
        "mesh2_heating_window_comparison": mesh_comparison,
        "next_action": "The selected 70:30 vol% target has now been compared at 1x and 2x through 0.2 ps; neither resolution melts. Preserve the old 41.6654 vol% results separately. Do not claim melt-pool convergence from this all-solid window; source property provenance and the experimental density convention remain open.",
        "property_limit": properties["source_scope"]["provenance_limit"],
    }
    (RUN_DIR / "composition_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
