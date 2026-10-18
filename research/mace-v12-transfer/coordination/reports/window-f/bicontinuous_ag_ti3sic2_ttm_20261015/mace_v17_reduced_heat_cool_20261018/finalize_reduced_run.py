#!/usr/bin/env python3
"""Offline aggregate/OVITO packaging for the reduced V17 heating run."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
MODEL_SHA = "559b12c024fe8e2155e26ffc7ca48cf9fa1c2682ea37ff9c5367af14a5cc43d3"
DT_FS = 0.5
MELT_CLUSTER_MIN = 20
MELT_LOWQ_FRACTION_MIN = 0.5
MELT_RMSD_MIN_A = 1.5


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def num(row, key):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def main():
    baseline = load_csv(ROOT / "BASELINE_300K_METRICS.csv")
    heating = load_csv(ROOT / "HEATING_METRICS.csv")
    geom = json.loads((ROOT / "REDUCED_GEOMETRY_AUDIT.json").read_text())
    static = json.loads((ROOT / "INITIAL_STATIC_REPORT.json").read_text())
    heat_rows = sorted(heating, key=lambda x: num(x, "time_ps"))
    base_rows = sorted(baseline, key=lambda x: num(x, "time_ps"))
    baseline_end_ps = num(base_rows[-1], "time_ps")
    # The segmented ASE restart code stored a cumulative heating clock, but
    # added it to the previous absolute time at each restart. Normalize from
    # the monotonic stage clock for analysis and OVITO metadata.
    for r in heat_rows:
        r["time_ps"] = str(baseline_end_ps + num(r, "stage_time_ps"))
    hot_base = [r for r in base_rows if int(float(r["n_hot_Ag"])) > 0]
    baseline_q6 = float(np.median([num(r, "hot_Ag_median_q6") for r in hot_base]))
    melt_hits = []
    for r in heat_rows:
        ok = (num(r, "T_hot_Ag_K") >= 1235.0 and
              num(r, "hot_Ag_RMSD_from_baseline_A") >= MELT_RMSD_MIN_A and
              num(r, "hot_Ag_median_q6") < 0.70 * baseline_q6 and
              num(r, "hot_Ag_low_q6_fraction") >= MELT_LOWQ_FRACTION_MIN and
              num(r, "hot_low_q6_largest_cluster_atoms") >= MELT_CLUSTER_MIN)
        melt_hits.append(ok)
    sustained_melt = any(all(melt_hits[i:i+3]) for i in range(max(0, len(melt_hits)-2)))
    max_cluster_row = max(heat_rows, key=lambda r: num(r, "hot_low_q6_largest_cluster_atoms"))
    max_temp_row = max(heat_rows, key=lambda r: num(r, "T_hot_Ag_K"))
    min_q_row = min(heat_rows, key=lambda r: num(r, "hot_Ag_median_q6"))
    max_rmsd_row = max(heat_rows, key=lambda r: num(r, "hot_Ag_RMSD_from_baseline_A"))
    max_force_row = max(heat_rows, key=lambda r: num(r, "Fmax_mobile_eV_A"))
    max_lowfrac_row = max(heat_rows, key=lambda r: num(r, "hot_Ag_low_q6_fraction"))
    end = heat_rows[-1]
    tsc_min_cc = min(num(r, "TSC_largest_component_fraction_cut3p25") for r in heat_rows)
    tsc_min_tic = min(num(r, "min_Ti_C_A") for r in heat_rows)

    # Preserve every recorded baseline/heating frame in one OVITO-ready extxyz.
    ovito = ROOT / "OVITO_HEAT_TRACE.extxyz"
    if ovito.exists():
        ovito.unlink()
    frame_count = 0
    for src in (ROOT / "BASELINE_300K.extxyz", ROOT / "HEATING.extxyz"):
        for atoms in read(src, index=":", format="extxyz"):
            if str(atoms.info.get("stage", "")).lower() == "heating":
                atoms.info["time_ps"] = baseline_end_ps + float(atoms.info.get("stage_time_ps", 0.0))
            write(ovito, atoms, format="extxyz", append=ovito.exists())
            frame_count += 1

    combined = base_rows + heat_rows
    if combined:
        fields = ["stage", "time_ps", "T_all_K", "T_Ag_K", "T_TSC_K", "T_hot_Ag_K", "T_cold_sink_K",
                  "Fmax_mobile_eV_A", "Fmax_Ag_eV_A", "hot_Ag_median_CN_3p45", "hot_Ag_median_q6",
                  "hot_Ag_low_q6_count", "hot_Ag_low_q6_fraction", "hot_low_q6_largest_cluster_atoms",
                  "hot_Ag_RMSD_from_baseline_A", "TSC_largest_component_fraction_cut3p25", "min_Ti_C_A",
                  "min_Ag_Ag_A", "min_Ag_TSC_A", "energy_eV"]
        with (ROOT / "TEMPERATURE_AND_STRUCTURE.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for r in combined:
                w.writerow(r)
        with (ROOT / "AG_LIQUID_CLUSTER_STATS.csv").open("w", newline="") as f:
            fields2 = ["stage", "time_ps", "T_hot_Ag_K", "hot_Ag_median_CN_3p45", "hot_Ag_median_q6",
                       "hot_Ag_low_q6_count", "hot_Ag_low_q6_fraction", "hot_low_q6_largest_cluster_atoms",
                       "hot_Ag_RMSD_from_baseline_A", "Fmax_Ag_eV_A", "TSC_largest_component_fraction_cut3p25",
                       "min_Ag_Ag_A", "min_Ag_TSC_A"]
            w = csv.DictWriter(f, fieldnames=fields2, extrasaction="ignore")
            w.writeheader()
            w.writerows(heat_rows)

    total_wall = (10.414138626001659 + 372.640006075002 + 1884.4039251020004 +
                  3803.6672978050046 + 3774.098505160997 + 9378.482604572004)
    result = {
        "MODEL_ATOMS": geom["atom_count"],
        "TIMESTEP_FS": DT_FS,
        "HEATING_PHYSICAL_TIME_PS": float(num(end, "stage_time_ps")),
        "COOLING_PHYSICAL_TIME_PS": 0.0,
        "AG_MELTING_ESTABLISHED": bool(sustained_melt),
        "AG_RESOLIDIFICATION_ESTABLISHED": False,
        "LIQUID_CLUSTER_MAX_ATOMS": int(num(max_cluster_row, "hot_low_q6_largest_cluster_atoms")),
        "LIQUID_CLUSTER_MAX_TIME_PS": num(max_cluster_row, "time_ps"),
        "MAX_HOT_Ag_LOW_Q6_FRACTION": num(max_lowfrac_row, "hot_Ag_low_q6_fraction"),
        "MAX_HOT_Ag_RMSD_A": num(max_rmsd_row, "hot_Ag_RMSD_from_baseline_A"),
        "MIN_HOT_Ag_Q6": num(min_q_row, "hot_Ag_median_q6"),
        "BASELINE_HOT_Ag_Q6": baseline_q6,
        "MAX_HOT_Ag_TEMPERATURE_K": num(max_temp_row, "T_hot_Ag_K"),
        "MAX_MOBILE_FORCE_eV_A": num(max_force_row, "Fmax_mobile_eV_A"),
        "TSC_BACKBONE_MIN_CONNECTED_FRACTION": tsc_min_cc,
        "TSC_MIN_Ti_C_DISTANCE_A": tsc_min_tic,
        "HEATING_STARTED": True,
        "COOLING_STARTED": False,
        "MELT_CRITERIA": {
            "three_consecutive_samples": True, "T_hot_Ag_min_K": 1235.0,
            "median_q6_below_fraction_of_baseline": 0.70, "low_q6_fraction_min": MELT_LOWQ_FRACTION_MIN,
            "RMSD_min_A": MELT_RMSD_MIN_A, "connected_low_q6_cluster_min_atoms": MELT_CLUSTER_MIN,
            "low_q6_threshold": 0.35,
        },
        "FINAL_HEATING_SAMPLE": {k: end[k] for k in ("time_ps", "T_hot_Ag_K", "hot_Ag_median_q6", "hot_Ag_low_q6_fraction",
                                                          "hot_low_q6_largest_cluster_atoms", "hot_Ag_RMSD_from_baseline_A",
                                                          "Fmax_mobile_eV_A", "TSC_largest_component_fraction_cut3p25")},
        "WALL_SECONDS_MEASURED_STATIC_BASELINE_HEATING": total_wall,
        "WALL_HOURS_MEASURED_STATIC_BASELINE_HEATING": total_wall / 3600.0,
        "trajectory_frames": frame_count,
        "OVITO_OUTPUT": str(ovito.name),
        "MACE_V17_sha256": MODEL_SHA,
        "geometry_mass_fraction_Ag": geom["Ag_mass_fraction"],
        "notes": [
            "The 1,084-atom system is a thin laterally periodic Ag(111)/Ti3SiC2(0001) interface slab, not the 39,486-atom bicontinuous target.",
            "First 0.30 ps used 300 mobile TSC atoms in the 300 K sink; the following 0.70 ps used a two-plane 200-atom TSC sink after the thin-substrate heat-transfer assessment.",
            "The Ag heat-source labels remained fixed at 94 atoms in the upper two Ag(111) planes within a 10 A lateral radius, target 1900 K.",
            "The peak connected low-q6 group never met the selected liquid-cluster, fraction, migration, and persistence criteria; heating was completed to 1.0 ps, but cooling was not started because no liquid Ag state was established.",
            "No cooling/resolidification claim is made. Provisional V17 and a small constrained slab make this exploratory evidence, not a validated physical melting point or production prediction.",
        ],
    }
    (ROOT / "THERMAL_RUN_STATUS.json").write_text(json.dumps(result, indent=2, allow_nan=True) + "\n")
    md = [
        "# Reduced MACE V17 local-heating run",
        "",
        "## Result",
        "",
        f"- Model atoms: {result['MODEL_ATOMS']} (Ti₃SiC₂={geom['TSC_atom_count']}, Ag={geom['Ag_atom_count']}); Ag mass fraction {100*geom['Ag_mass_fraction']:.2f}%. ",
        f"- Time step: {DT_FS:.1f} fs; local-heating duration: {result['HEATING_PHYSICAL_TIME_PS']:.3f} ps; cooling duration: 0 ps.",
        f"- Heating wall time including static force and 300 K baseline: {total_wall/3600:.2f} h ({total_wall:.1f} s).",
        f"- Ag melt established: {'YES' if sustained_melt else 'NO'}; maximum connected low-q₆ Ag cluster: {result['LIQUID_CLUSTER_MAX_ATOMS']} atoms at {result['LIQUID_CLUSTER_MAX_TIME_PS']:.3f} ps.",
        f"- Maximum low-q₆ fraction: {100*result['MAX_HOT_Ag_LOW_Q6_FRACTION']:.1f}%; minimum hot-region median q₆ {result['MIN_HOT_Ag_Q6']:.3f} (300 K baseline {baseline_q6:.3f}); maximum hot-region RMS displacement {result['MAX_HOT_Ag_RMSD_A']:.3f} Å.",
        f"- Maximum mobile force {result['MAX_MOBILE_FORCE_eV_A']:.2f} eV/Å; minimum TSC connectivity fraction {tsc_min_cc:.3f}.",
        "",
        "The trajectory reached 1 ps under a 1900 K local Langevin target, but the measured Ag remained predominantly ordered: q₆ did not fall below 70% of its baseline for three samples, fewer than half of hot Ag atoms had q₆ < 0.35, the largest connected low-q₆ cluster was below 20 atoms, and RMS migration stayed below 1.5 Å. Therefore this run shows heating and local disorder only, not a demonstrated melt pool. Cooling was deliberately not run because no liquid state was established.",
        "",
        "## Boundary and model limits",
        "",
        "The reduced geometry keeps the original 30.575 × 26.479 Å lateral coincidence cell and ideal Ag/Ti₃SiC₂ interface. It contains a 7.54 Å six-plane stoichiometric Ti₃SiC₂ slab, four Ag(111) planes, 18 Å top vacuum, XY periodic boundaries, only the bottom Ti plane fixed, and all other TSC atoms mobile. It is a layered interface patch rather than the original bicontinuous 39,486-atom model. Ag/TSC mass balance shifts to 72.73 wt% Ag.",
        "",
        "The first 0.30 ps used a 300-atom TSC heat sink. To reduce thermostat clamping in this thin slab, the remaining 0.70 ps used only the two mobile TSC planes immediately above the fixed support (200 atoms) at 300 K; other atoms remained dynamically mobile. The same V17 weights and 0.5 fs step were used throughout.",
        "",
        "## Files",
        "",
        "- `OVITO_HEAT_TRACE.extxyz`: baseline plus 1 ps heating trajectory; no cooling frames.",
        "- `TEMPERATURE_AND_STRUCTURE.csv`: time series of temperature, force, q₆, coordination, displacement proxy, contacts, and TSC connectivity.",
        "- `AG_LIQUID_CLUSTER_STATS.csv`: per-frame Ag disorder and connected-cluster indicators.",
        "- `HEATING.extxyz`, `HEATING_METRICS.csv`, and checkpoint files: raw staged heating outputs.",
        "- `REDUCED_MODEL.extxyz` / `REDUCED_MODEL.data`: independent reduced slab coordinates.",
        "",
        "This result is a MACE V17 exploratory trajectory. No DFT, retraining, potential mixing, or MD with EAM/TSCBOP/Morse was used. The potential remains provisional; non-melting at this finite time does not establish an equilibrium melting temperature.",
        "",
    ]
    (ROOT / "FINAL_EXPLORATORY_REPORT.md").write_text("\n".join(md))

    # Hash all deliverables except the checksum manifest itself.
    files = sorted(p for p in ROOT.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt")
    with (ROOT / "SHA256SUMS.txt").open("w") as f:
        for path in files:
            f.write(f"{sha(path)}  {path.name}\n")
    print(json.dumps(result, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
