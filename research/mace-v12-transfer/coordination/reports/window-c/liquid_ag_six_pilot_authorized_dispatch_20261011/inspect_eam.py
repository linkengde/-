#!/usr/bin/env python3
"""Read an EAM parameter file and metadata; never calculate energy or forces.

An unavailable or unverified file is not replaced with a different potential.
ASE 3.29 parses funcfl (eam), setfl (alloy), and Finnis--Sinclair (fs).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys


FORM_ALIASES = {"setfl": "alloy", "funcfl": "eam", "alloy": "alloy",
                "eam": "eam", "fs": "fs"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def detect_format(path: Path, requested: str = "auto") -> str:
    if requested != "auto":
        if requested not in FORM_ALIASES:
            raise ValueError("UNSUPPORTED_EAM_FORMAT; use eam, alloy, fs, setfl or funcfl")
        return FORM_ALIASES[requested]
    suffix = path.suffix.lower()
    if suffix in {".alloy", ".setfl"}:
        return "alloy"
    if suffix == ".fs":
        return "fs"
    if suffix == ".eam":
        # Some published setfl files use .eam: inspect the header, not the name.
        lines = path.read_text(encoding="utf-8").splitlines()
        if len(lines) >= 5:
            fields = lines[3].split()
            try:
                count = int(fields[0])
            except (ValueError, IndexError):
                count = -1
            if count > 0 and len(fields) == count + 1:
                return "alloy"
        return "eam"
    raise ValueError("FORMAT_NOT_FROZEN; specify the published file's actual format")


def inspect_potential(path: Path, form: str, box_edges: list[float] | None = None):
    """Parse splines only; caller is responsible for source/license release."""
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("POTENTIAL_NOT_AVAILABLE")
    import numpy as np
    from ase.calculators.eam import EAM
    import ase

    if ase.__version__ != "3.29.0":
        raise ValueError(f"ASE_VERSION_MISMATCH: expected 3.29.0, got {ase.__version__}")
    actual_form = detect_format(path, form)
    calc = EAM(potential=str(path), form=actual_form)
    if "Ag" not in calc.elements:
        raise ValueError("POTENTIAL_HAS_NO_AG")
    arrays = [calc.embedded_data, calc.density_data, calc.rphi_data,
              calc.mass, calc.Z, calc.r, calc.rho]
    if not all(np.isfinite(array).all() for array in arrays):
        raise ValueError("POTENTIAL_CONTAINS_NAN_OR_INF")
    if (calc.nr < 4 or calc.nrho < 4 or calc.dr <= 0 or calc.drho <= 0
            or not math.isfinite(calc.cutoff) or calc.cutoff <= 0):
        raise ValueError("INVALID_EAM_GRIDS_OR_CUTOFF")
    ag_index = calc.elements.index("Ag")
    if int(calc.Z[ag_index]) != 47 or not 107.0 < float(calc.mass[ag_index]) < 109.0:
        raise ValueError("AG_ELEMENT_METADATA_MISMATCH")
    last_r = float(calc.r[-1])
    gap = max(0.0, float(calc.cutoff) - last_r)
    tail_radii = np.linspace(last_r, float(calc.cutoff), 11) if gap else np.asarray([float(calc.cutoff)])
    density_function = (calc.electron_density[ag_index, ag_index] if actual_form == "fs"
                        else calc.electron_density[ag_index])
    density_derivative = (calc.d_electron_density[ag_index, ag_index] if actual_form == "fs"
                          else calc.d_electron_density[ag_index])
    pair_function = calc.phi[ag_index, ag_index]
    pair_derivative = calc.d_phi[ag_index, ag_index]
    # Table/interpolant metadata only: no Atoms or total E/F evaluation.
    density_tail = density_function(tail_radii)
    phi_tail = pair_function(tail_radii)
    density_table = (calc.density_data[ag_index, ag_index] if actual_form == "fs"
                     else calc.density_data[ag_index])
    tail = {
        "grid_last_r_A": last_r, "nr_times_dr_A": float(calc.nr * calc.dr),
        "gap_A": gap, "requires_explicit_C_policy": bool(gap > 1e-12),
        "within_one_dr": bool(gap <= calc.dr * (1 + 1e-8) + 1e-12),
        "sample_count": len(tail_radii), "sample_r_A": tail_radii.tolist(),
        "max_abs_density": float(np.max(np.abs(density_tail))),
        "max_abs_phi_eV": float(np.max(np.abs(phi_tail))),
        "density_at_sample_r": np.asarray(density_tail).tolist(),
        "phi_eV_at_sample_r": np.asarray(phi_tail).tolist(),
        "density_derivative_at_sample_r": np.asarray(density_derivative(tail_radii)).tolist(),
        "phi_derivative_eV_per_A_at_sample_r": np.asarray(pair_derivative(tail_radii)).tolist(),
        "last_five_table_r_A": calc.r[-5:].tolist(),
        "last_five_density_table_values": density_table[-5:].tolist(),
        "last_five_rphi_table_eV_A": calc.rphi_data[ag_index, ag_index, -5:].tolist(),
        "status": "C_TAIL_REVIEW_REQUIRED" if gap > 1e-12 else "NO_RADIAL_EXTRAPOLATION_REQUIRED",
        "scope": "PARAMETER_INTERPOLANT_METADATA_NO_ATOMIC_ENERGY_FORCE_OR_MD",
    }
    if not all(np.isfinite(values).all() for values in (
            density_tail, phi_tail, density_derivative(tail_radii), pair_derivative(tail_radii))):
        raise ValueError("PARAMETER_RADIAL_TAIL_NOT_FINITE")
    boxes = []
    for edge in box_edges or []:
        if not math.isfinite(edge) or edge <= 0:
            raise ValueError("INVALID_BOX_EDGE")
        boxes.append({"box_edge_A": edge, "half_box_A": edge / 2,
                      "cutoff_exceeds_half_box": bool(calc.cutoff > edge / 2),
                      "multiple_periodic_images_required": bool(calc.cutoff > edge / 2),
                      "RDF_max_radius_A": edge / 2})
    metadata = {
        "operation": "READ_PARAMETER_METADATA_ONLY_NO_ENERGY_FORCE_OR_MD",
        "file_name": path.name, "file_bytes": path.stat().st_size,
        "potential_sha256": sha256_file(path), "ase_version": ase.__version__,
        "format": actual_form, "elements_in_file_order": list(calc.elements),
        "atomic_numbers": np.asarray(calc.Z).astype(int).tolist(),
        "masses_amu": np.asarray(calc.mass).tolist(), "Ag_element_index": ag_index,
        "cutoff_A": float(calc.cutoff), "nr": int(calc.nr), "dr_A": float(calc.dr),
        "radial_grid_max_A": float(calc.r[-1]), "nrho": int(calc.nrho),
        "radial_tail": tail,
        "drho": float(calc.drho), "density_grid_max": float(calc.rho[-1]),
        "header_lines": [line.rstrip("\n") for line in calc.header],
        "units": "LAMMPS metal/EAM table convention: eV, Angstrom, amu; source must confirm",
        "neighbor_implementation": "ASE PrimitiveNeighborList with integer offsets; bothways=True",
        "minimum_image_only": False, "periodic_image_support": "SOURCE_REVIEW_PASS_NOT_NUMERICAL_RUN",
        "box_checks": boxes, "license_verified_by_this_parser": False,
        "source_verified_by_this_parser": False,
    }
    return metadata, calc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--potential", type=Path, required=True)
    parser.add_argument("--format", choices=["auto", *FORM_ALIASES], default="auto")
    parser.add_argument("--box-edge-A", type=float, action="append", default=[])
    parser.add_argument("--expected-sha256")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        if not args.potential.is_file():
            raise ValueError("POTENTIAL_NOT_AVAILABLE; no replacement or download is attempted")
        if args.expected_sha256 and sha256_file(args.potential) != args.expected_sha256.lower():
            raise ValueError("POTENTIAL_SHA_MISMATCH")
        metadata, _ = inspect_potential(args.potential, args.format, args.box_edge_A)
        body = json.dumps(metadata, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
        if args.output:
            if args.output.exists():
                raise ValueError("OUTPUT_EXISTS; preserve the previous inspection")
            args.output.write_text(body, encoding="utf-8")
        print(body, end="")
        return 0
    except Exception as exc:
        print(json.dumps({"status": "INSPECTION_HOLD", "reason": str(exc)},
                         ensure_ascii=False), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
