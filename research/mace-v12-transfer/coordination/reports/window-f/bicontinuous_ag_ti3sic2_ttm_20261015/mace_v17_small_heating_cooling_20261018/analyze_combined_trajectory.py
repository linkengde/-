#!/usr/bin/env python3
"""Offline multi-indicator Ag disorder, migration, cooldown and TSC audit."""
from pathlib import Path
import csv, json, sys
import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import analyze_ag_structure as ag

baseline = read(ROOT / "BASELINE_300K.extxyz", index=":", format="extxyz")
reference = baseline[0]
ref_ids, _, ref_q6, _ = ag.ag_order(reference)
base_records = []
ag.summarize_frames(baseline, "baseline", reference, ref_q6, base_records)

sources = [("heating_1fs", "HEAT_COOL.extxyz"),
           ("heating_0p2fs", "HEATING_DT0P2FS.extxyz"),
           ("cooling_0p2fs", "COOLING_DT0P2FS.extxyz")]
records = []
for name, fn in sources:
    path = ROOT / fn
    if path.exists():
        frames = read(path, index=":", format="extxyz")
        ag.summarize_frames(frames, name, reference, ref_q6, records)

hot_ref = [r for r in base_records if r["hot_Ag_atoms"]]
baseline_q6 = float(np.median([r["hot_Ag_median_q6"] for r in hot_ref])) if hot_ref else float("nan")
heat_rows = [r for r in records if r["trajectory"].startswith("heating")]
cool_rows = [r for r in records if r["trajectory"].startswith("cooling")]
hottest = max(heat_rows, key=lambda r: r["T_hot_Ag_K"] if np.isfinite(r["T_hot_Ag_K"]) else -1) if heat_rows else None
coolend = cool_rows[-1] if cool_rows else None
melt = bool(hottest and hottest["T_hot_Ag_K"] >= 1235 and
            hottest["hot_Ag_RMSD_from_300K_start_A"] >= 1.5 and
            hottest["hot_Ag_median_q6"] < 0.70 * baseline_q6 and
            hottest["hot_Ag_low_q6_fraction"] >= 0.5 and
            hottest["largest_hot_low_q6_Ag_component_atoms"] >= 20)
reorder = bool(melt and coolend and
               coolend["T_hot_Ag_K"] < hottest["T_hot_Ag_K"] - 100 and
               coolend["hot_Ag_median_q6"] > hottest["hot_Ag_median_q6"] + 0.10 and
               coolend["hot_Ag_low_q6_fraction"] < hottest["hot_Ag_low_q6_fraction"])
all_rows = base_records + records
with (ROOT / "AG_MELTING_COOLING_ANALYSIS.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(all_rows[0]))
    w.writeheader(); w.writerows(all_rows)
report = {
    "q6_cutoff_A": ag.CUTOFF,
    "low_q6_tag_threshold": 0.35,
    "hot_region_definition": "upper three Ag(111) planes within 8 A lateral radius of cell center",
    "baseline_hot_Ag_median_q6": baseline_q6,
    "hottest_heating_frame": hottest,
    "final_cooling_frame": coolend,
    "liquid_candidate_test": "hot Ag T>=1235 K, hot-region RMSD>=1.5 A, q6 median <70% of baseline, >=50% q6<0.35, largest connected low-q6 Ag cluster >=20 atoms",
    "AG_MELTING_EVIDENCE": "POSSIBLE_LOCAL_MELT_CANDIDATE" if melt else "NOT_ESTABLISHED",
    "AG_RESOLIDIFICATION_EVIDENCE": "POSSIBLE_REORDERING_CANDIDATE" if reorder else "NOT_ESTABLISHED",
    "TSC_backbone_min_connectivity_fraction": min((r["TSC_largest_component_fraction_cut3p25"] for r in records), default=float("nan")),
    "maximum_framewise_min_Ti_C_A": max((r["TSC_min_Ti_C_A"] for r in records), default=float("nan")),
    "interpretation": "Thermal vibration, elevated temperature, or an isolated close contact do not establish melting. A cooling segment without preceding melt evidence does not establish resolidification.",
}
(ROOT / "AG_MELTING_COOLING_ANALYSIS.json").write_text(json.dumps(report, indent=2, allow_nan=True) + "\n")
print(json.dumps(report, allow_nan=True))
