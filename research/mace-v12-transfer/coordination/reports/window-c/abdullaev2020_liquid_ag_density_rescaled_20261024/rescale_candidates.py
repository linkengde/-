#!/usr/bin/env python3
"""Create separate geometry-only Ag snapshots anchored to Abdullaev et al. (2020)."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read, write


HERE = Path(__file__).resolve().parent
PREP = HERE.parent / "continuation_geometry_prep_20261024"
OUT = HERE / "candidates"
MASS_G_MOL = 107.8682
AVOGADRO = 6.02214076e23
SOURCE_PDF_SHA256 = "cd37e2b91be0369ebe1b9d72ffb7377a28dfc74171ed0dd3cd22892f95e23c3a"
SOURCE_DOI = "10.1088/1742-6596/1677/1/012161"
SOURCE_EQUATION = "rho_m(T) = 9327.1 - 0.9956 * (T - 1234.93) kg m^-3"
OLD_ANCHOR = 8.867447227420516
RATIOS = (0.97, 1.00, 1.03)
TEMPERATURES = (1250, 1600)


def density_g_cm3(natoms: int, volume_a3: float) -> float:
    return natoms * MASS_G_MOL / AVOGADRO * 1.0e24 / volume_a3


def liquid_reference_g_cm3(temperature_k: int) -> float:
    if not 1234.93 <= temperature_k <= 1631:
        raise ValueError(f"Equation (2) is outside its stated liquid range at {temperature_k} K")
    return (9327.1 - 0.9956 * (temperature_k - 1234.93)) / 1000.0


def source_paths() -> list[Path]:
    paths = list((PREP / "ag_liquid_pilot").glob("Ag_liquid_*/*.extxyz"))
    paths += list((PREP / "ag_liquid_extended_candidates").glob("Ag_liquid_*/*.extxyz"))
    return sorted(
        path for path in paths
        if any(f"Ag_liquid_{temperature}K_" in path.name for temperature in TEMPERATURES)
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    records: list[dict] = []
    for source_path in source_paths():
        source = read(source_path)
        info = dict(source.info)
        temperature = int(info["temperature_K"])
        ratio = float(info["density_ratio_to_model_anchor"])
        source_config = str(info["config_type"])
        if temperature not in TEMPERATURES or ratio not in RATIOS:
            raise ValueError(f"Unexpected candidate metadata: {source_path}")

        source_density = density_g_cm3(len(source), source.get_volume())
        old_expected = OLD_ANCHOR * ratio
        if abs(source_density - old_expected) > 1.0e-5:
            raise ValueError(f"Source density metadata mismatch in {source_path}")

        reference_density = liquid_reference_g_cm3(temperature)
        target_density = reference_density * ratio
        linear_scale = (source_density / target_density) ** (1.0 / 3.0)
        scaled = source.copy()
        fractional_before = scaled.get_scaled_positions(wrap=False)
        scaled.set_cell(np.asarray(scaled.cell) * linear_scale, scale_atoms=True)
        fractional_after = scaled.get_scaled_positions(wrap=False)
        fractional_delta = (fractional_after - fractional_before + 0.5) % 1.0 - 0.5
        max_fractional_change = float(np.max(np.abs(fractional_delta)))
        actual_density = density_g_cm3(len(scaled), scaled.get_volume())
        if abs(actual_density - target_density) > 1.0e-8:
            raise ValueError(f"Rescaled density mismatch in {source_path}")
        if max_fractional_change > 1.0e-10:
            raise ValueError(f"Fractional coordinates changed in {source_path}")

        distances = scaled.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        minimum_distance = float(np.min(distances))
        if (
            len(scaled) != 64
            or scaled.get_chemical_formula() != "Ag64"
            or not np.isfinite(scaled.positions).all()
            or not np.isfinite(np.asarray(scaled.cell)).all()
            or scaled.get_volume() <= 0
            or minimum_distance <= 1.8
        ):
            raise ValueError(f"Invalid geometry after scaling in {source_path}")

        old_model_anchor = float(info.pop("model_density_anchor_g_cm3"))
        old_ratio = float(info.pop("density_ratio_to_model_anchor"))
        config = f"{source_config}_Abdullaev2020ref"
        info.update({
            "config_type": config,
            "candidate_role": "B1_LIQ_training_candidate_pending_role_audit",
            "dataset_role": "geometry_only_proposal_no_DFT_label",
            "proposal_only_no_DFT_label": True,
            "dft_launch_enabled": False,
            "density_reference_source": "Abdullaev, Khairulin & Stankus (2020), J. Phys.: Conf. Ser. 1677, 012161",
            "density_reference_doi": SOURCE_DOI,
            "density_reference_equation": SOURCE_EQUATION,
            "density_reference_valid_range_K": "1234.93-1631",
            "density_reference_temperature_K": temperature,
            "density_reference_value_g_cm3": reference_density,
            "density_ratio_to_literature_reference": ratio,
            "density_target_g_cm3": target_density,
            "legacy_model_density_anchor_g_cm3": old_model_anchor,
            "legacy_density_ratio_to_model_anchor": old_ratio,
            "legacy_source_configuration": source_config,
            "isotropic_linear_scale_from_legacy_geometry": linear_scale,
            "source_pdf_sha256": SOURCE_PDF_SHA256,
            "density_reference_uncertainty_note": (
                "1250 K: neighboring Table 1 density errors are 0.25% at Tf and 0.26% at 1300 K; "
                "coverage/confidence level is unspecified. 1600 K: Table 1 reports 0.34%; "
                "coverage/confidence level is unspecified. Density-ratio offsets are a sensitivity grid, not uncertainty."
            ),
            "geometry_transformation": "isotropic cell/coordinate scaling; fractional coordinates preserved",
        })
        scaled.info.clear()
        scaled.info.update(info)

        output_dir = OUT / config
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{config}.extxyz"
        write(output_path, scaled, format="extxyz")
        written = read(output_path)
        written_density = density_g_cm3(len(written), written.get_volume())
        written_fractional = (written.get_scaled_positions(wrap=False) - fractional_before + 0.5) % 1.0 - 0.5
        if (
            written.get_chemical_formula() != "Ag64"
            or len(written) != 64
            or abs(written_density - target_density) > 1.0e-8
            or np.max(np.abs(written_fractional)) > 2.0e-9
            or written.info["proposal_only_no_DFT_label"] is not True
            or written.info["dft_launch_enabled"] is not False
        ):
            raise ValueError(f"Serialized candidate failed read-back checks: {output_path}")
        records.append({
            "config_type": config,
            "source_file": str(source_path.relative_to(PREP)),
            "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "file": str(output_path.relative_to(HERE)),
            "file_sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
            "atoms": len(scaled),
            "formula": scaled.get_chemical_formula(),
            "temperature_K": temperature,
            "reference_density_g_cm3": reference_density,
            "density_ratio": ratio,
            "target_density_g_cm3": target_density,
            "actual_density_g_cm3": actual_density,
            "serialized_density_g_cm3": written_density,
            "linear_scale_factor": linear_scale,
            "max_fractional_coordinate_change": max_fractional_change,
            "minimum_pair_distance_A": minimum_distance,
            "no_DFT_label": scaled.info["proposal_only_no_DFT_label"],
            "DFT_launch_enabled": scaled.info["dft_launch_enabled"],
        })

    if len(records) != 12 or len({r["file_sha256"] for r in records}) != 12:
        raise ValueError(f"Expected 12 unique candidates, got {len(records)}")
    manifest = {
        "status": "PASS_GEOMETRY_ONLY_RESCALED_CANDIDATES",
        "candidate_count": len(records),
        "source": {
            "citation": "R. N. Abdullaev, R. A. Khairulin, and S. V. Stankus, J. Phys.: Conf. Ser. 1677 (2020) 012161",
            "doi": SOURCE_DOI,
            "license": "Creative Commons Attribution 3.0 (as stated in the article PDF)",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "equation_2": SOURCE_EQUATION,
            "liquid_validity_range_K": [1234.93, 1631],
        },
        "method": "Isotropically rescaled copies; fractional coordinates preserved; source candidates were not modified.",
        "legacy_model_density_anchor_g_cm3": OLD_ANCHOR,
        "density_ratios": list(RATIOS),
        "temperatures_included_K": list(TEMPERATURES),
        "temperature_excluded": {
            "1050 K": "Below the 1234.93 K melting point; the paper's liquid equation is out of range, so no liquid-reference rescaling was performed."
        },
        "limitations": [
            "Geometry-only proposals; no DFT labels, no training integration, and no phase/property validation claim.",
            "The source formula supplies fitted reference densities, not measurements at every tabulated temperature.",
            "At 1250 K the paper's adjacent tabulated density-error entries are 0.25% and 0.26%; no confidence level is specified.",
            "The 0.97/1.00/1.03 density ratios are sensitivity proposals and are not measurement uncertainty bounds.",
        ],
        "records": records,
    }
    manifest_path = HERE / "rescaled_candidate_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "status": manifest["status"],
        "candidate_count": len(records),
        "density_range_g_cm3": [min(r["actual_density_g_cm3"] for r in records), max(r["actual_density_g_cm3"] for r in records)],
        "minimum_pair_distance_A": min(r["minimum_pair_distance_A"] for r in records),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
