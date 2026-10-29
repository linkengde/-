#!/usr/bin/env python3
"""Run low-cost TTM parameter-limit checks for the archived 1 ps fields.

This performs reduced ODE sensitivity calculations only; it does not launch
or modify any production simulation.
"""
from __future__ import annotations
import csv, hashlib, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
MODEL = REPO / "models/ag_ti3sic2_qsgs_1to1_finer"
PHASE = MODEL / "phase_field_qsgs"
DATA = HERE / "source_data"
sys.path.insert(0, str(PHASE))
from qsgs_source_integration import integrate_qsgs_source_power_W

EV_J = 1.602176634e-19
A3_M3 = 1.0e-30


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def weighted_quantile(x: np.ndarray, w: np.ndarray, q: float) -> float:
    order = np.argsort(x)
    return float(x[order[min(np.searchsorted(np.cumsum(w[order]), q), len(x)-1)]])


def integrate_exchange(q: np.ndarray, weights: np.ndarray, ttab: np.ndarray,
                       gtab: np.ndarray, uej: np.ndarray, cpv: float,
                       gscale: float, dt_s: float, nsteps: int) -> dict:
    ue = np.zeros_like(q, dtype=float)
    hl = np.zeros_like(q, dtype=float)
    te = np.full_like(q, 300.0, dtype=float)
    tl = np.full_like(q, 300.0, dtype=float)
    for _ in range(nsteps):
        g = np.interp(te, ttab, gtab) * gscale
        exchange = g * (te - tl)
        ue = np.maximum(ue + (q - exchange) * dt_s, 0.0)
        hl += exchange * dt_s
        te = np.interp(ue, uej, ttab)
        tl = 300.0 + hl / cpv
    source_j_m3 = float(np.dot(weights, q) * dt_s * nsteps)
    total_j_m3 = float(np.dot(weights, ue + hl))
    return {
        "Te_mean_K": float(np.dot(weights, te)),
        "Te_q99_K": weighted_quantile(te, weights, .99),
        "Te_q999_K": weighted_quantile(te, weights, .999),
        "Tl_mean_K": float(np.dot(weights, tl)),
        "Tl_q99_K": weighted_quantile(tl, weights, .99),
        "Tl_q999_K": weighted_quantile(tl, weights, .999),
        "lattice_energy_fraction": float(np.dot(weights, hl) / source_j_m3),
        "energy_residual_relative": (total_j_m3-source_j_m3)/source_j_m3,
        "max_Te_K": float(te.max()),
        "source_energy_density_J_m3": source_j_m3,
    }


def main() -> None:
    geom_path = MODEL / "qsgs_morphology.npz"
    table_path = MODEL / "ttm_data/Ag_TTM_hybrid_figure_fit_300_20000K.csv"
    summary_path = DATA / "qsgs_component_native_heating_to_1ps.json"
    progress_path = DATA / "qsgs_component_native_heating_to_1ps_recovery_progress.json"
    fields_path = DATA / "qsgs_component_native_heating_to_1ps_fields.npz"
    with np.load(geom_path, allow_pickle=False) as z:
        base_ag = ~z["skeleton"].astype(bool)
        size0 = np.asarray(z["size_A"], dtype=float)
        seed = int(z["seed"])
    ag = np.tile(base_ag, (3, 3, 1))
    size_A = size0 * np.array([3., 3., 1.])
    spacing_A = size_A / np.asarray(ag.shape)
    voxel_A3 = float(np.prod(spacing_A))
    voxel_m3 = voxel_A3 * A3_M3
    rows = list(csv.DictReader(table_path.open(encoding="utf-8-sig")))
    ttab = np.asarray([float(r["Te_K"]) for r in rows])
    cetab = np.asarray([float(r["Ce_LinFig2c_J_m-3_K-1"]) for r in rows])
    gtab = np.asarray([float(r["G_LinFig2d_W_m-3_K-1"]) for r in rows])
    ketab = np.asarray([float(r["ke_KarnaFig4b_W_m-1_K-1"]) for r in rows])
    uej = np.asarray([float(r["Ue_delta_from_300K_eV_A-3"]) for r in rows]) * EV_J / A3_M3
    with np.load(fields_path, allow_pickle=False) as z:
        te = z["electron_temperature_K"]
        tl = z["lattice_temperature_K"]
    ag_te_mean = float(te[ag].mean())
    ag_tl_mean = float(tl[ag].mean())
    rho = 10490.0
    cp_mass = 25.4 / 0.1078682
    cp_volume = rho * cp_mass

    # Same 1 ps integrated source, uniformly mixed over the Ag phase.
    peak_flux = 5.0 * EV_J / (1.0e-12 * 1.0e-20)
    source_cell_W = integrate_qsgs_source_power_W(ag, ag.shape, size_A, peak_flux, 30.0, 20.0)
    p_total = float(source_cell_W.sum())
    ag_volume_m3 = float(ag.sum() * voxel_m3)
    q_uniform = p_total / ag_volume_m3
    dt_s, nsteps = 1.0e-16, 10_000
    hom = []
    for scale in (.5, 1., 2.):
        result = integrate_exchange(np.asarray([q_uniform]), np.asarray([1.]), ttab, gtab, uej,
                                    cp_volume, scale, dt_s, nsteps)
        hom.append({"G_multiplier": scale, **result,
                    "model": "uniformly mixed 0D; same 1 ps total source; no spatial gradients"})

    # Extreme no-electron-diffusion limit; source remains spatially localized.
    q_cells = np.sort(source_cell_W[ag] / voxel_m3)
    nbin = min(4096, q_cells.size)
    groups = np.array_split(q_cells, nbin)
    q_bins = np.asarray([group.mean() for group in groups])
    weights = np.asarray([group.size for group in groups], dtype=float)
    weights /= weights.sum()
    isolated = []
    for scale in (.5, 1., 2.):
        result = integrate_exchange(q_bins, weights, ttab, gtab, uej, cp_volume,
                                    scale, dt_s, nsteps)
        isolated.append({"G_multiplier": scale, **result,
                         "source_rate_bins": nbin,
                         "model": "independent Ag voxels, electron diffusion off; 4096 equal-count bins; reduced diagnostic, not full graph rerun"})

    ce_at_mean = float(np.interp(ag_te_mean, ttab, cetab))
    ke_at_mean = float(np.interp(ag_te_mean, ttab, ketab))
    g_at_mean = float(np.interp(ag_te_mean, ttab, gtab))
    alpha = ke_at_mean / ce_at_mean
    times_s = np.asarray([.1, .5, 1.]) * 1e-12
    ke_sensitivity = []
    for factor in (.5, 1., 2.):
        a = factor * alpha
        ke_sensitivity.append({
            "ke_multiplier": factor,
            "alpha_m2_s": a,
            "sqrt_alpha_t_nm_at_0p1_0p5_1ps": (np.sqrt(a*times_s)*1e9).tolist(),
            "note": "analytic scaling; not a full 3D ke-perturbed run",
        })

    summary = json.loads(summary_path.read_text())
    comparison = {
        "audited_1ps_geometry": {
            "path": "models/ag_ti3sic2_qsgs_1to1_finer/qsgs_morphology.npz",
            "sha256": sha256(geom_path), "seed": seed,
            "shape_xyz": list(base_ag.shape), "Ag_volume_fraction": float(base_ag.mean()),
            "Ti3SiC2_volume_fraction": float(1-base_ag.mean()),
            "estimated_Ag_mass_fraction_using_10p49_and_4p52_g_cm3": float(
                base_ag.mean()*10.49 / (base_ag.mean()*10.49 + (1-base_ag.mean())*4.52)
            ),
            "estimated_Ag_to_Ti3SiC2_mass_ratio_using_those_densities": float(
                base_ag.mean()*10.49 / ((1-base_ag.mean())*4.52)
            ),
        },
        "later_distinct_geometry_not_used_in_this_audit": {
            "path": "models/ag_ti3sic2_qsgs_1to1_finer/runs/single_temperature_mesh2_validation_20261027/ag607wt_ti3sic2_393wt_nano_20261028/nano_geometry_ag607wt_ti3sic2_393wt.npz",
            "sha256": "edebee1d8d254a60e51c7de0b63d96c10263020989c9c963449b932323e95695",
            "Ag_volume_fraction": .399586, "Ti3SiC2_volume_fraction": .600414,
            "basis": "60.7/39.3 wt% converted with candidate densities; thresholded existing QSGS growth-rank realization, not a regenerated independent realization",
            "parent_geometry_sha256": "7e0e1509b414597ce4e2a1c9fc387ed9656caf507f1186e2d28de79c3c1bfc43",
        },
    }
    report = {
        "baseline": {
            "source_power_W_recomputed": p_total,
            "source_energy_1ps_J": p_total*1e-12,
            "observed_source_energy_J": summary["samples"][-1]["source_energy_J"],
            "observed_Te_mean_K": ag_te_mean, "observed_Tl_mean_K": ag_tl_mean,
            "observed_liquid_volume_fraction": summary["samples"][-1]["ag_liquid_volume_fraction"],
            "energy_residual_relative": summary["samples"][-1]["energy_residual_relative_to_source"],
            "geometry": comparison["audited_1ps_geometry"],
            "geometry_distinction": comparison["later_distinct_geometry_not_used_in_this_audit"],
        },
        "thermal_parameter_scales_at_observed_mean_Te": {
            "Te_mean_K": ag_te_mean, "Ce_J_m3_K": ce_at_mean,
            "ke_W_m_K": ke_at_mean, "G_W_m3_K": g_at_mean,
            "alpha_m2_s": alpha, "Ce_over_G_ps": ce_at_mean/g_at_mean*1e12,
            "rhoCp_over_G_ps": cp_volume/g_at_mean*1e12,
            "lattice_rhoCp_over_G_note": "lattice response time for a fixed 1-K Te-Tl gap; not the coupled equilibration time",
            "property_range_K": [float(ttab.min()), float(ttab.max())],
        },
        "ke_scale_analytic_sensitivity": ke_sensitivity,
        "G_scale_homogeneous_0D_sensitivity": hom,
        "no_electron_diffusion_localized_source_sensitivity": isolated,
        "method_notes": [
            "0D tests use tabulated Ue(T) inversion and G(T), constant room-temperature Ag lattice volumetric heat capacity, dt=0.1 fs for 1 ps; no spatial conduction.",
            "No-diffusion test uses exact cell-integrated source power on the Ag mask and independent-cell response grouped into 4096 equal-count source-rate bins.",
            "ke scaling is sqrt(ke) scaling of the one-sided characteristic diffusion length sqrt(alpha*t), not a 3D parameter-perturbed simulation.",
            "These are controlled reduced diagnostics, not predictions that replace the archived full graph run.",
        ],
        "hashes": {str(p.relative_to(REPO)): sha256(p) for p in [geom_path, table_path, summary_path, progress_path, fields_path]},
    }
    out_json = HERE / "sensitivity_results.json"
    out_json.write_text(json.dumps(report, indent=2)+"\n")
    out_csv = HERE / "sensitivity_results.csv"
    with out_csv.open("w", newline="") as f:
        fields = ["case", "multiplier", "Te_mean_K", "Te_q99_K", "Te_q999_K", "Tl_mean_K", "Tl_q99_K", "Tl_q999_K", "lattice_energy_fraction", "energy_residual_relative", "max_Te_K", "note"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for r in hom:
            w.writerow({"case":"well_mixed_0D_G", "multiplier":r["G_multiplier"], "Te_mean_K":r["Te_mean_K"], "Tl_mean_K":r["Tl_mean_K"], "lattice_energy_fraction":r["lattice_energy_fraction"], "energy_residual_relative":r["energy_residual_relative"], "max_Te_K":r["max_Te_K"], "note":r["model"]})
        for r in isolated:
            w.writerow({"case":"no_diffusion_local_source_G", "multiplier":r["G_multiplier"], "Te_mean_K":r["Te_mean_K"], "Te_q99_K":r["Te_q99_K"], "Te_q999_K":r["Te_q999_K"], "Tl_mean_K":r["Tl_mean_K"], "Tl_q99_K":r["Tl_q99_K"], "Tl_q999_K":r["Tl_q999_K"], "lattice_energy_fraction":r["lattice_energy_fraction"], "energy_residual_relative":r["energy_residual_relative"], "max_Te_K":r["max_Te_K"], "note":r["model"]})
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
