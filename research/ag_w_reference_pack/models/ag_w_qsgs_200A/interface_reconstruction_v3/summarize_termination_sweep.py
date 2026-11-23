#!/usr/bin/env python3
"""Aggregate all preregistered v3 termination/registry runs into one audit report."""
from __future__ import annotations

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RUNS = HERE / "termination_runs"
CELL_MATCH = json.loads((HERE / "cell_match_report.json").read_text())
DFT = {
    "Ag100_W100": {"work_of_separation_J_m2": 2.80, "interface_energy_J_m2": 0.72},
    "Ag111_W110": {"work_of_separation_J_m2": 2.45, "interface_energy_J_m2": 0.86},
}
DFT_TOL = 0.30


def main() -> None:
    result = {
        "title": "Ag-W v3 low-strain one-interface termination and registry sweep",
        "protocol": "V3_PROTOCOL.md",
        "cell_match_report": "cell_match_report.json",
        "thresholds_changed_after_calculation": False,
        "dft_screening_tolerance_J_m2": DFT_TOL,
        "exploratory_use_scope": "A literature mismatch does not by itself forbid qualitative graduate-level simulation. The local static mechanics gate is reported separately; quantitative interface energetics are not validated.",
        "orientations": {},
    }
    overall_pass = True
    for orientation in DFT:
        reports = []
        for path in sorted(RUNS.glob(f"{orientation}_tAg*_tW*.json")):
            report = json.loads(path.read_text())
            if "potential_sha256" not in result:
                result["potential_sha256"] = report.get("potential_sha256")
                result["lammps_version"] = report.get("lammps_version")
            item = report["orientations"][orientation]
            best = item.get("best_interface_summary")
            registries = item.get("registries", {})
            valid_summaries = [r["summary"] for r in registries.values() if "summary" in r]
            scan_points = [point for reg in registries.values() for point in reg["scan"]]
            reports.append(
                {
                    "file": path.name,
                    "termination_indices": item["termination_indices"],
                    "best_registry": item.get("best_registry"),
                    "best_summary": best,
                    "registry_summaries": {
                        name: row["summary"] for name, row in registries.items() if "summary" in row
                    },
                    "run_count": len(scan_points),
                    "successful_run_count": sum(bool(point.get("run_ok")) for point in scan_points),
                    "force_converged_count": sum(
                        bool(point.get("run_ok")) and point.get("stopping_criterion") == "force tolerance"
                        for point in scan_points
                    ),
                    "maximum_mobile_force_eV_A": max(
                        (float(point["max_mobile_force_eV_A"]) for point in scan_points if point.get("run_ok")),
                        default=None,
                    ),
                    "minimum_distances_over_all_scans_A": {
                        pair: min(
                            float(point["minimum_distances_A"][pair])
                            for point in scan_points
                            if point.get("run_ok") and point.get("minimum_distances_A", {}).get(pair) is not None
                        )
                        for pair in ("Ag_Ag", "W_W", "Ag_W")
                    },
                    "same_phase_collapse_count_over_all_scans": sum(
                        point.get("run_ok")
                        and (
                            point.get("minimum_distances_A", {}).get("Ag_Ag", 1.0e9)
                            < 0.85 * (4.0900219316 / (2.0 ** 0.5))
                            or point.get("minimum_distances_A", {}).get("W_W", 1.0e9)
                            < 0.85 * (3.0 ** 0.5 * 3.1599967891 / 2.0)
                        )
                        for point in scan_points
                    ),
                    "atom_count_mismatch_count": sum(
                        point.get("run_ok") and point.get("retained_atom_count") != point.get("atom_count")
                        for point in scan_points
                    ),
                    "registry_minimum_energy_range_eV": (
                        max(s["minimum_energy_eV"] for s in valid_summaries)
                        - min(s["minimum_energy_eV"] for s in valid_summaries)
                        if valid_summaries
                        else None
                    ),
                }
            )
        if not reports:
            raise RuntimeError(f"No termination runs found for {orientation}")
        expected_terms = 4 if orientation == "Ag100_W100" else 6
        complete_terms = len(reports) == expected_terms
        all_points = [r for config in reports for r in config["registry_summaries"].values()]
        global_best_config = min(reports, key=lambda r: r["best_summary"]["minimum_energy_eV"])
        best = global_best_config["best_summary"]
        target = DFT[orientation]
        dft_pass = (
            abs(best["work_of_separation_error_vs_DFT_J_m2"]) <= DFT_TOL
            and abs(best["interface_energy_error_vs_DFT_J_m2"]) <= DFT_TOL
            and best["interface_energy_J_m2"] > 0
        )
        mechanical_pass = (
            best["force_gate_pass"]
            and not best["same_phase_collapse_flag"]
            and best["atom_count_preserved"]
        )
        item = {
            "selected_match": CELL_MATCH["orientations"][orientation]["selected_cell"],
            "dft_reference": target,
            "termination_combinations_expected": expected_terms,
            "termination_combinations_completed": len(reports),
            "termination_coverage_complete": complete_terms,
            "registry_translation_count_per_termination": 3 if orientation == "Ag100_W100" else 4,
            "gap_points_per_registry": 10,
            "total_interface_minimizations": sum(r["run_count"] for r in reports),
            "successful_interface_minimizations": sum(r["successful_run_count"] for r in reports),
            "force_converged_interface_minimizations": sum(r["force_converged_count"] for r in reports),
            "maximum_mobile_force_over_all_scans_eV_A": max(r["maximum_mobile_force_eV_A"] for r in reports),
            "minimum_distances_over_all_scans_A": {
                pair: min(r["minimum_distances_over_all_scans_A"][pair] for r in reports)
                for pair in ("Ag_Ag", "W_W", "Ag_W")
            },
            "same_phase_collapse_count_over_all_scans": sum(
                r["same_phase_collapse_count_over_all_scans"] for r in reports
            ),
            "atom_count_mismatch_count": sum(r["atom_count_mismatch_count"] for r in reports),
            "global_lowest_energy_configuration": {
                "termination_indices": global_best_config["termination_indices"],
                "registry": global_best_config["best_registry"],
                "summary": best,
            },
            "registry_energy_spread_by_termination_eV": {
                f"Ag{r['termination_indices']['Ag']}_W{r['termination_indices']['W']}": r[
                    "registry_minimum_energy_range_eV"
                ]
                for r in reports
            },
            "mechanical_gate_pass": mechanical_pass,
            "DFT_screening_gate_pass": dft_pass,
            "orientation_gate_pass": complete_terms and mechanical_pass and dft_pass,
            "exploratory_static_use_pass": complete_terms and mechanical_pass,
            "termination_results": reports,
        }
        result["orientations"][orientation] = item
        overall_pass = overall_pass and item["orientation_gate_pass"]
    result["status"] = "LOCAL_INTERFACE_SCREEN_PASSED" if overall_pass else "LOCAL_INTERFACE_SCREEN_REJECTED"
    result["exploratory_static_use_assessment"] = (
        "MECHANICALLY_SUFFICIENT_FOR_QUALITATIVE_LOCAL_INTERFACE_EXPLORATION"
        if all(row["exploratory_static_use_pass"] for row in result["orientations"].values())
        else "NOT_MECHANICALLY_READY"
    )
    (HERE / "local_interface_v3_report.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "status": result["status"],
        "orientations": {
            name: {
                "termination_coverage_complete": row["termination_coverage_complete"],
                "interface_minimizations": row["total_interface_minimizations"],
                "force_converged": row["force_converged_interface_minimizations"],
                "global_best": row["global_lowest_energy_configuration"],
                "DFT_screening_gate_pass": row["DFT_screening_gate_pass"],
            }
            for name, row in result["orientations"].items()
        },
    }, indent=2))


if __name__ == "__main__":
    main()
