#!/usr/bin/env python3
"""Join the V10 5 A scan with the V9 10 A post-cold reference scan."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
MODEL=HERE.parent
scan5=json.loads((HERE/"v10_5A_phase_scan.json").read_text())
scan10=json.loads((MODEL/"interface_reconstruction_v9_hybrid"/"v9_postcold_10A_phase_scan.json").read_text())
mapping=json.loads((HERE/"mapping_5A_best_final.json").read_text())
transfer=json.loads((HERE/"v10_transfer_distance_sensitivity.json").read_text())
best5=min(scan5["phase_results"],key=lambda x:x["production_capacity_flow"]["locally_unassignable_overflow_A3"])
best10=scan10["best_by_shared_capacity_flow"]
phase10=next(x for x in scan10["phase_results"] if x["offset_xy_A"]==best10["offset_xy_A"])
profile10=transfer["raw_source_profile"]
result={
    "title":"Ag-W V10 5 A versus retained 10 A phase comparison",
    "frozen_v9_commit":scan5["frozen_v9_commit"],
    "input_restart_sha256":scan5["input"]["restart_sha256"],
    "comparison":[
        {"grid_A":5.0,"phase_count":scan5["phases_tested"],
         "best_phase_xy_A":best5["offset_xy_A"],
         "raw_max_fill":best5["raw_fill_max_combined"],
         "raw_overflow_demand_A3":best5["raw_overflow_demand_A3"],
         "locally_unassignable_A3_at_20A":best5["production_capacity_flow"]["locally_unassignable_overflow_A3"],
         "source_power_eV_ps":best5["source_spatial_profile"]["power_eV_ps"],
         "source_radial_RMS_A":best5["source_spatial_profile"]["radial_RMS_A"],
         "source_mean_depth_A":best5["source_spatial_profile"]["mean_depth_A"],
         "all_phases_pass":scan5["all_phases_pass"],
         "real_lammps_mapping_only_status":mapping["status"],
         "real_lammps_unassignable_A3":best5["production_capacity_flow"]["locally_unassignable_overflow_A3"]},
        {"grid_A":10.0,"phase_count":scan10["phase_count"],
         "best_phase_xy_A":best10["offset_xy_A"],
         "raw_max_fill":phase10["raw_fill_max_combined"],
         "raw_overflow_demand_A3":best10["overflow_demand_A3"],
         "locally_unassignable_A3_at_20A":best10["unassigned_overflow_A3"],
         "source_power_eV_ps":profile10["power_eV_ps"],
         "source_radial_RMS_A":profile10["radial_RMS_A"],
         "source_mean_depth_A":profile10["mean_depth_A"],
         "all_phases_pass":False,
         "real_lammps_mapping_only_status":"NOT_RUN_IN_V10; V9 previously failed mapping preflight"}
    ],
    "interpretation":"5 A decreases the best-phase residual relative to the 10 A V9 scan, but all four 5 A phases fail and the actual best-phase LAMMPS mapping-only preflight fails. Residuals depend on cell volume/resolution and are not a pass criterion by themselves.",
    "no_ttm_or_source_step":True}
(HERE/"v10_resolution_comparison.json").write_text(json.dumps(result,indent=2)+"\n")
