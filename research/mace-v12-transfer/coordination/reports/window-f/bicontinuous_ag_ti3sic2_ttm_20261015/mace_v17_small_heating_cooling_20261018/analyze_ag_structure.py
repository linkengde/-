#!/usr/bin/env python3
"""Offline Ag structure/order and mobility summaries for the V17 trajectory."""
from pathlib import Path
import csv
import json
import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list
from scipy.special import sph_harm_y
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parent
CUTOFF = 3.45
L = 6


def ag_order(atoms):
    symbols = np.asarray(atoms.get_chemical_symbols())
    ids = np.flatnonzero(symbols == "Ag")
    i, j, D = neighbor_list("ijD", atoms, CUTOFF)
    ag = (symbols[i] == "Ag") & (symbols[j] == "Ag")
    i, j, D = i[ag], j[ag], D[ag]
    deg = np.bincount(i, minlength=len(atoms))
    qlm = np.zeros((len(atoms), 2 * L + 1), dtype=complex)
    if len(i):
        r = np.linalg.norm(D, axis=1)
        theta = np.arccos(np.clip(D[:, 2] / r, -1, 1))
        phi = np.arctan2(D[:, 1], D[:, 0])
        for k, m in enumerate(range(-L, L + 1)):
            y = sph_harm_y(L, m, theta, phi)
            qlm[:, k] = np.bincount(i, weights=y.real, minlength=len(atoms)) + 1j * np.bincount(i, weights=y.imag, minlength=len(atoms))
        valid = deg > 0
        qlm[valid] /= deg[valid, None]
    q6 = np.zeros(len(atoms))
    q6[ids] = np.sqrt(4 * np.pi / (2 * L + 1) * np.sum(np.abs(qlm[ids]) ** 2, axis=1))
    cn = deg[ids]
    # A connected component of low-q6 Ag gives a spatial-contiguity indicator,
    # not an independent thermodynamic definition of liquid.
    disordered = ids[q6[ids] < 0.35]
    index = np.full(len(atoms), -1, dtype=int)
    index[disordered] = np.arange(len(disordered))
    ee = (index[i] >= 0) & (index[j] >= 0)
    adj = coo_matrix((np.ones(2 * np.sum(ee)),
                      (np.r_[index[i[ee]], index[j[ee]]], np.r_[index[j[ee]], index[i[ee]]])),
                     shape=(len(disordered), len(disordered))).tocsr()
    if len(disordered):
        _, labels = connected_components(adj, directed=False)
        sizes = np.bincount(labels)
        largest = int(sizes.max())
    else:
        largest = 0
    return ids, cn, q6[ids], largest


def largest_component(atoms, selected_ids):
    if len(selected_ids) == 0:
        return 0
    symbols = np.asarray(atoms.get_chemical_symbols())
    i, j = neighbor_list("ij", atoms, CUTOFF)
    allowed = np.zeros(len(atoms), dtype=bool)
    allowed[selected_ids] = True
    keep = allowed[i] & allowed[j] & (symbols[i] == "Ag") & (symbols[j] == "Ag")
    index = np.full(len(atoms), -1, dtype=int)
    index[selected_ids] = np.arange(len(selected_ids))
    adj = coo_matrix((np.ones(2 * np.sum(keep)),
                      (np.r_[index[i[keep]], index[j[keep]]], np.r_[index[j[keep]], index[i[keep]]])),
                     shape=(len(selected_ids), len(selected_ids))).tocsr()
    _, labels = connected_components(adj, directed=False)
    return int(np.bincount(labels).max())


def backbone_metrics(atoms):
    phase = np.asarray(atoms.arrays["phase"])
    symbols = np.asarray(atoms.get_chemical_symbols())
    ids = np.flatnonzero(phase == "Ti3SiC2")
    index = np.full(len(atoms), -1, dtype=int)
    index[ids] = np.arange(len(ids))
    i, j, d = neighbor_list("ijd", atoms, 3.25)
    keep = (index[i] >= 0) & (index[j] >= 0)
    adj = coo_matrix((np.ones(2 * np.sum(keep)),
                      (np.r_[index[i[keep]], index[j[keep]]],
                       np.r_[index[j[keep]], index[i[keep]]])),
                     shape=(len(ids), len(ids))).tocsr()
    _, labels = connected_components(adj, directed=False)
    tsc_cc_fraction = float(np.bincount(labels).max() / len(ids)) if len(ids) else 0.0
    undirected = i < j
    tic = undirected & (((symbols[i] == "Ti") & (symbols[j] == "C")) |
                        ((symbols[i] == "C") & (symbols[j] == "Ti")))
    return {
        "TSC_largest_component_fraction_cut3p25": tsc_cc_fraction,
        "TSC_min_Ti_C_A": float(np.min(d[tic])) if np.any(tic) else float("nan"),
    }


def lateral_delta(atoms, reference):
    cell = np.asarray(atoms.cell.array)
    frac = atoms.positions @ np.linalg.inv(cell)
    frac0 = reference.positions @ np.linalg.inv(np.asarray(reference.cell.array))
    duv = frac[:, :2] - frac0[:, :2]
    dz = (frac[:, 2] - frac0[:, 2]) * cell[2, 2]
    best = np.full((len(atoms), 2), np.inf)
    bestnorm = np.full(len(atoms), np.inf)
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            xy = (duv + np.array([a, b])) @ cell[:2, :2]
            norm = np.linalg.norm(xy, axis=1)
            take = norm < bestnorm
            best[take] = xy[take]
            bestnorm[take] = norm[take]
    return np.column_stack([best, dz])


def hot_mask(atoms, ids):
    layer = np.asarray(atoms.arrays["Ag111_layer"], dtype=int)
    frac = atoms.positions @ np.linalg.inv(np.asarray(atoms.cell.array))
    duv = frac[:, :2] - np.array([0.5, 0.5])
    best = np.full(len(atoms), np.inf)
    cell = np.asarray(atoms.cell.array)
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            xy = (duv + np.array([a, b])) @ cell[:2, :2]
            best = np.minimum(best, np.linalg.norm(xy, axis=1))
    return (layer[ids] >= 3) & (best[ids] <= 8.0)


def summarize_frames(frames, stage_name, reference, baseline_q6, records):
    refids, _, _, _ = ag_order(reference)
    ref_hot = hot_mask(reference, refids)
    # Keep the diagnostic region tied to the same initial Ag atoms. Re-selecting
    # by a moving spatial radius can silently change the atom count during heat
    # and cooling and confound order/mobility comparisons.
    fixed_hot_ids = refids[ref_hot]
    ref_order = float(np.median(baseline_q6[ref_hot]))
    refpos = reference.copy()
    ref_symbols = np.asarray(reference.get_chemical_symbols())
    for frame in frames:
        ids, cn, q6, largest = ag_order(frame)
        is_hot = np.isin(ids, fixed_hot_ids)
        hot_low_q6_ids = ids[is_hot & (q6 < 0.35)]
        hot_low_q6_component = largest_component(frame, hot_low_q6_ids)
        disp = lateral_delta(frame, refpos)[ids]
        if np.any(is_hot):
            displacement = disp[is_hot]
            displacement -= displacement.mean(axis=0)
            rmsd = float(np.sqrt(np.mean(np.sum(displacement ** 2, axis=1))))
        else:
            rmsd = float("nan")
        phase = np.asarray(frame.arrays["phase"])
        support = np.asarray(frame.arrays["support_fixed"], dtype=bool)
        force_norm = np.asarray(frame.arrays["mace_force_norm_eV_A"], dtype=float)
        atom_ids = np.asarray(frame.arrays["id"], dtype=int)
        global_max_index = int(np.argmax(force_norm))
        hot_force_norm = force_norm[ids[is_hot]]
        tsc_ids = np.flatnonzero((phase == "Ti3SiC2") & ~support)
        tsc_disp = lateral_delta(frame, refpos)[tsc_ids]
        if len(tsc_ids):
            tsc_disp -= tsc_disp.mean(axis=0)
            tsc_rmsd = float(np.sqrt(np.mean(np.sum(tsc_disp ** 2, axis=1))))
            tsc_dmax = float(np.max(np.linalg.norm(tsc_disp, axis=1)))
        else:
            tsc_rmsd = tsc_dmax = float("nan")
        tsc_geom = backbone_metrics(frame)
        hot_indices = ids[is_hot]
        hot_velocities = frame.get_velocities()[hot_indices]
        hot_masses = frame.get_masses()[hot_indices, None]
        hot_vcm = np.sum(hot_masses * hot_velocities, axis=0) / np.sum(hot_masses)
        hot_dof = max(1, 3 * len(hot_indices) - 3)
        hot_temperature = float(np.sum(hot_masses * (hot_velocities - hot_vcm) ** 2) /
                                (hot_dof * 8.617333262145e-5))
        low_hot_ids = ids[is_hot & (q6 < 0.35)]
        if len(low_hot_ids):
            pos = frame.positions[low_hot_ids]
            spans = np.ptp(pos, axis=0)
        else:
            spans = np.array([0.0, 0.0, 0.0])
        info = frame.info
        records.append({
            "trajectory": stage_name,
            "stage": str(info.get("stage", "unknown")),
            "step": int(info.get("step", -1)),
            "time_ps": float(info.get("time_ps", 0.0)),
            "hot_Ag_atoms": int(np.sum(is_hot)),
            "T_hot_Ag_K": hot_temperature,
            "T_Ag_K": float(info.get("T_Ag_K", float("nan"))),
            "Fmax_all_eV_A": float(force_norm[global_max_index]),
            "Fmax_all_atom_id": int(atom_ids[global_max_index]),
            "Fmax_all_element": str(frame.get_chemical_symbols()[global_max_index]),
            "Fmax_all_phase": str(phase[global_max_index]),
            "Fmax_hot_Ag_eV_A": float(np.max(hot_force_norm)) if len(hot_force_norm) else float("nan"),
            "hot_Ag_median_CN_3p45A": float(np.median(cn[is_hot])) if np.any(is_hot) else float("nan"),
            "hot_Ag_median_q6": float(np.median(q6[is_hot])) if np.any(is_hot) else float("nan"),
            "hot_Ag_low_q6_fraction": float(np.mean(q6[is_hot] < 0.35)) if np.any(is_hot) else float("nan"),
            "hot_Ag_low_q6_count": int(np.sum(q6[is_hot] < 0.35)),
            "hot_Ag_q6_vs_300K_baseline": float(np.median(q6[is_hot]) - ref_order) if np.any(is_hot) else float("nan"),
            "hot_Ag_RMSD_from_300K_start_A": rmsd,
            "largest_hot_low_q6_Ag_component_atoms": hot_low_q6_component,
            "hot_low_q6_region_span_x_A": float(spans[0]),
            "hot_low_q6_region_span_y_A": float(spans[1]),
            "hot_low_q6_region_span_z_A": float(spans[2]),
            "TSC_mobile_RMSD_from_300K_start_A": tsc_rmsd,
            "TSC_mobile_max_displacement_A": tsc_dmax,
            **tsc_geom,
        })


def main():
    baseline = read(ROOT / "BASELINE_300K.extxyz", index=":", format="extxyz")
    ref = baseline[0]
    rids, _, rq6, _ = ag_order(ref)
    baseline_rows = []
    records = []
    summarize_frames(baseline, "baseline", ref, rq6, baseline_rows)
    trajectory_path = ROOT / "HEAT_COOL.extxyz"
    if trajectory_path.exists():
        heatcool = read(trajectory_path, index=":", format="extxyz")
        summarize_frames(heatcool, "heat_cool", ref, rq6, records)
    all_rows = baseline_rows + records
    if all_rows:
        with (ROOT / "AG_STRUCTURE_CLASSIFICATION.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0]))
            w.writeheader(); w.writerows(all_rows)
    base_hot = [r for r in baseline_rows if r["hot_Ag_atoms"]]
    heat = [r for r in records if r["stage"] == "heating"]
    cool = [r for r in records if r["stage"] == "cooling"]
    base_q = float(np.median([r["hot_Ag_median_q6"] for r in base_hot])) if base_hot else float("nan")
    hottest = max(heat, key=lambda r: r["T_hot_Ag_K"] if np.isfinite(r["T_hot_Ag_K"]) else -1) if heat else None
    coolest_end = cool[-1] if cool else None
    melt_candidate = bool(hottest and hottest["T_hot_Ag_K"] >= 1235 and
                          hottest["hot_Ag_RMSD_from_300K_start_A"] >= 1.5 and
                          hottest["hot_Ag_median_q6"] < 0.70 * base_q and
                          hottest["hot_Ag_low_q6_fraction"] >= 0.5 and
                          hottest["largest_hot_low_q6_Ag_component_atoms"] >= 20)
    resolid_candidate = bool(melt_candidate and coolest_end and
                             coolest_end["T_hot_Ag_K"] < (hottest["T_hot_Ag_K"] - 100) and
                             coolest_end["hot_Ag_median_q6"] > hottest["hot_Ag_median_q6"] + 0.10 and
                             coolest_end["hot_Ag_low_q6_fraction"] < hottest["hot_Ag_low_q6_fraction"])
    report = {
        "neighbor_cutoff_A": CUTOFF,
        "order_parameter": "per-atom Steinhardt q6 magnitude from Ag-Ag neighbors; q6<0.35 tagged locally disordered for cluster bookkeeping",
        "reference": "first frame of the 300 K baseline trajectory; RMSD is relative to this frame with lateral PBC and hot-region COM removed",
        "baseline_hot_Ag_median_q6": base_q,
        "highest_temperature_heating_frame": hottest,
        "last_cooling_frame": coolest_end,
        "AG_MELTING_EVIDENCE": "POSSIBLE_LOCAL_MELT_CANDIDATE" if melt_candidate else "NOT_ESTABLISHED",
        "AG_RESOLIDIFICATION_EVIDENCE": "POSSIBLE_REORDERING_CANDIDATE" if resolid_candidate else "NOT_ESTABLISHED",
        "criterion_note": "These conservative multi-indicator flags are exploratory screens only; sub-ps trajectories cannot by themselves establish equilibrium melting or a reliable welded joint.",
    }
    (ROOT / "AG_STRUCTURE_CLASSIFICATION_REPORT.json").write_text(json.dumps(report, indent=2, allow_nan=True) + "\n")
    print(json.dumps(report, allow_nan=True))


if __name__ == "__main__":
    main()
