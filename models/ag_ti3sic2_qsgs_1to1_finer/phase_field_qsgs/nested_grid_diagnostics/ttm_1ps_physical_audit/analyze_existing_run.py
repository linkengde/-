#!/usr/bin/env python3
"""Post-process an existing native-QSGS 1 ps TTM checkpoint; does not run MD/TTM."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

EV_J = 1.602176634e-19
A3_M3 = 1.0e-30


def weighted_stats(a: np.ndarray, mask: np.ndarray) -> dict:
    x = np.asarray(a[mask], dtype=float)
    q = np.quantile(x, [0.0, 0.1, 0.5, 0.9, 0.99, 0.999, 1.0])
    return {
        "mean_K": float(x.mean()), "std_K": float(x.std()),
        "q0_K": float(q[0]), "q10_K": float(q[1]), "q50_K": float(q[2]),
        "q90_K": float(q[3]), "q99_K": float(q[4]), "q999_K": float(q[5]),
        "max_K": float(q[6]),
    }


def circular_extent(values: np.ndarray, n: int, spacing: float) -> tuple[int, int, float]:
    u = np.unique(values)
    if u.size == 0:
        return 0, -1, 0.0
    if u.size == n:
        return 0, n - 1, n * spacing
    gaps = (np.roll(u, -1) - u) % n
    empty = gaps - 1
    gi = int(np.argmax(empty))
    start = int(u[(gi + 1) % u.size])
    span = int(n - empty[gi])
    return start, (start + span - 1) % n, span * spacing


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--model-dir", type=Path, required=True)
    p.add_argument("--diag-dir", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    geom_path = args.model_dir / "qsgs_morphology.npz"
    fields_path = args.diag_dir / "qsgs_component_native_heating_to_1ps_fields.npz"
    run_path = args.diag_dir / "qsgs_component_native_heating_to_1ps.json"
    progress_path = args.diag_dir / "qsgs_component_native_heating_to_1ps_recovery_progress.json"
    table_path = args.model_dir / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv"
    with np.load(geom_path, allow_pickle=False) as g:
        base_ag = ~np.asarray(g["skeleton"], dtype=bool)
        base_size_A = np.asarray(g["size_A"], dtype=float)
        seed = int(g["seed"])
    ag = np.tile(base_ag, (3, 3, 1))
    size_A = base_size_A * np.array([3.0, 3.0, 1.0])
    spacing_A = size_A / np.asarray(ag.shape)
    voxel_A3 = float(np.prod(spacing_A))
    with np.load(fields_path, allow_pickle=False) as s:
        te = np.asarray(s["electron_temperature_K"], dtype=float)
        tl = np.asarray(s["lattice_temperature_K"], dtype=float)
        fl = np.asarray(s["ag_liquid_fraction"], dtype=float)
        ue = np.asarray(s["electron_energy_density_eV_A3"], dtype=float)
        hl = np.asarray(s["lattice_enthalpy_density_J_m3"], dtype=float)
        cum_cpl = np.asarray(s["electron_to_lattice_cumulative_energy_density_J_m3"], dtype=float)
    if te.shape != ag.shape:
        raise ValueError(f"geometry/field shape mismatch: {ag.shape} vs {te.shape}")
    run = json.loads(run_path.read_text())
    progress = json.loads(progress_path.read_text())
    final = run["samples"][-1]
    times = [0.1, 0.2, 0.5, 1.0]
    samples = []
    for t in times:
        matches = [r for r in progress["samples"] if abs(float(r["time_ps"]) - t) <= 1e-8]
        if matches:
            r = matches[-1]
            samples.append({"time_ps": t, "Te_mean_K": r["electron_temperature_mean_K"],
                            "Te_max_K": r["electron_temperature_max_K"],
                            "Tl_mean_K": r["lattice_temperature_mean_K"],
                            "Tl_max_K": r["lattice_temperature_max_K"],
                            "liquid_volume_fraction": r["ag_liquid_volume_fraction"],
                            "energy_residual_relative": r["energy_residual_relative_to_source"]})
    with (args.out_dir / "temperature_history_0p1_0p2_0p5_1ps.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(samples[0]))
        w.writeheader(); w.writerows(samples)

    te_stats = weighted_stats(te, ag)
    tl_stats = weighted_stats(tl, ag)
    idx = np.unravel_index(np.argmax(np.where(ag, tl, -np.inf)), tl.shape)
    max_xyz_A = [(idx[k] + 0.5) * spacing_A[k] for k in range(3)]
    source_xy = size_A[:2] / 2.0
    dx = (max_xyz_A[0] - source_xy[0] + size_A[0] / 2) % size_A[0] - size_A[0] / 2
    dy = (max_xyz_A[1] - source_xy[1] + size_A[1] / 2) % size_A[1] - size_A[1] / 2
    depth_A = size_A[2] - max_xyz_A[2]

    def zone_metrics(field: np.ndarray, frac: float) -> dict:
        peak_excess = max(float(field[ag].max()) - 300.0, 0.0)
        active = ag & (field >= 300.0 + frac * peak_excess)
        coords = np.argwhere(active)
        if not coords.size:
            return {"threshold_fraction_of_peak_excess": frac, "voxel_count": 0, "volume_A3": 0.0}
        x0, x1, lx = circular_extent(coords[:, 0], ag.shape[0], spacing_A[0])
        y0, y1, ly = circular_extent(coords[:, 1], ag.shape[1], spacing_A[1])
        z0, z1 = int(coords[:, 2].min()), int(coords[:, 2].max())
        return {
            "threshold_fraction_of_peak_excess": frac,
            "threshold_K": 300.0 + frac * peak_excess,
            "voxel_count": int(coords.shape[0]),
            "volume_A3": float(coords.shape[0] * voxel_A3),
            "volume_fraction_of_Ag": float(coords.shape[0] / ag.sum()),
            "extent_A_xyz": [lx, ly, (z1-z0+1)*spacing_A[2]],
            "voxel_bbox_xyz": [[x0,x1],[y0,y1],[z0,z1]],
        }

    # Reconstruct per-voxel integrated source power and cumulative local
    # electron diffusion from the saved energy fields: dE_e = Q - E_e->l + D_e.
    from qsgs_source_integration import integrate_qsgs_source_power_W
    peak_flux_W_m2 = 5.0 * EV_J / (1.0e-12 * 1.0e-20)
    source_W = integrate_qsgs_source_power_W(ag, ag.shape, size_A, peak_flux_W_m2, 30.0, 20.0)
    duration_s = 1.0e-12
    voxel_volume_m3 = voxel_A3 * A3_M3
    electron_J_cell = ue * EV_J * voxel_A3
    source_J_cell = source_W * duration_s
    coupling_J_cell = cum_cpl * voxel_volume_m3
    diffusion_J_cell = electron_J_cell - source_J_cell + coupling_J_cell
    diffusion_net_J = float(diffusion_J_cell[ag].sum())
    diffusion_gross_J = float(np.abs(diffusion_J_cell[ag]).sum() / 2.0)

    # Instantaneous final-time graph-equivalent electron heat currents across
    # actual shared Ag voxel faces. x/y faces are periodic; z faces are internal
    # only. This is reconstructed from the final field, not archived flux history.
    rows = list(csv.DictReader(table_path.open(encoding="utf-8-sig")))
    temp_table = np.array([float(r["Te_K"]) for r in rows])
    ke_table = np.array([float(r["ke_KarnaFig4b_W_m-1_K-1"]) for r in rows])
    ke = np.interp(te, temp_table, ke_table)
    hot_te = ag & (te >= 300.0 + 0.5 * (float(te[ag].max()) - 300.0))
    face_flux_abs_W = 0.0
    hot_net_out_W = 0.0
    axis_gross_W = []
    for axis in range(3):
        Ti = te
        Tj = np.roll(te, -1, axis=axis)
        mi = ag
        mj = np.roll(ag, -1, axis=axis)
        ki = ke
        kj = np.roll(ke, -1, axis=axis)
        face = mi & mj
        if axis == 2:
            face[:, :, -1] = False
        kh = np.zeros_like(ke)
        kh[face] = 2 * ki[face] * kj[face] / (ki[face] + kj[face])
        area = np.prod(np.delete(spacing_A, axis)) * 1.0e-20
        dist = spacing_A[axis] * 1.0e-10
        flux_i_to_j = kh * (area / dist) * (Ti - Tj)
        flux_i_to_j[~face] = 0.0
        gross = float(np.abs(flux_i_to_j).sum())
        out = float(np.sum(flux_i_to_j * (hot_te.astype(float) - np.roll(hot_te, -1, axis=axis).astype(float))))
        axis_gross_W.append(gross)
        face_flux_abs_W += gross
        hot_net_out_W += out

    cp_300 = 25.4 / 0.1078682
    rho = 10490.0
    cl = rho * cp_300
    temp_at_mean = float(te_stats["mean_K"])
    Ce = float(np.interp(temp_at_mean, temp_table, [float(r["Ce_LinFig2c_J_m-3_K-1"]) for r in rows]))
    G = float(np.interp(temp_at_mean, temp_table, [float(r["G_LinFig2d_W_m-3_K-1"]) for r in rows]))
    K = float(np.interp(temp_at_mean, temp_table, ke_table))
    alpha = K / Ce
    timescale_s = np.array([0.1, 0.5, 1.0]) * 1.0e-12
    diffusion_lengths_nm = np.sqrt(alpha * timescale_s) * 1e9
    electron_tau_ps = Ce / G * 1e12
    lattice_tau_ps = cl / G * 1e12

    cp_molar=25.4
    ag_mass = ag.sum() * voxel_A3 * A3_M3 * rho
    source_total = float(source_W.sum() * duration_s)
    predicted_mean_dT_all_to_lattice = source_total / (ag_mass * cp_300)
    predicted_mean_dT_coupled = float(final["lattice_enthalpy_J"]) / (ag_mass * cp_300)

    # Source neighborhood: within one lateral Gaussian radius and one skin depth.
    x_center, y_center = size_A[:2] / 2.0
    x_pos = (np.arange(ag.shape[0]) + 0.5) * spacing_A[0]
    y_pos = (np.arange(ag.shape[1]) + 0.5) * spacing_A[1]
    z_pos = (np.arange(ag.shape[2]) + 0.5) * spacing_A[2]
    dx_src = np.minimum(np.abs(x_pos - x_center), size_A[0] - np.abs(x_pos - x_center))[:, None, None]
    dy_src = np.minimum(np.abs(y_pos - y_center), size_A[1] - np.abs(y_pos - y_center))[None, :, None]
    depth_src = (size_A[2] - z_pos)[None, None, :]
    source_neighborhood = ag & (dx_src**2 + dy_src**2 <= 30.0**2) & (depth_src <= 20.0)

    metrics = {
        "provenance": {
            "geometry_file": "models/ag_ti3sic2_qsgs_1to1_finer/qsgs_morphology.npz", "geometry_seed": seed,
            "geometry_source_shape_xyz": list(base_ag.shape),
            "geometry_ag_volume_fraction": float(base_ag.mean()),
            "tiled_domain_xyz_voxels": list(ag.shape),
            "tiled_domain_size_A": size_A.tolist(),
            "voxel_spacing_A": spacing_A.tolist(),
            "ag_native_voxels": int(ag.sum()), "ag_total_volume_A3": float(ag.sum()*voxel_A3),
            "ag_mass_kg": float(ag_mass), "tile_xy": 3,
            "summary": "source_data/" + run_path.name, "field_snapshot": "source_data/" + fields_path.name,
            "boundary_conditions": run["boundary_conditions"],
            "I0_eV_ps_A2": run["heating"]["peak_flux_eV_ps_A2"],
            "R_A": run["heating"]["gaussian_radius_A"], "gaussian": run["heating"]["gaussian_definition"],
            "skin_depth_A": run["heating"]["absorption_depth_A"],
            "actual_integrated_source_power_W_from_mask": float(source_W.sum()),
        },
        "1ps_temperature_statistics_ag_weighted": {"Te": te_stats, "Tl": tl_stats},
        "lattice_max_location": {
            "voxel_index_xyz": [int(v) for v in idx], "position_A_center": max_xyz_A,
            "distance_xy_from_source_center_A": float(np.hypot(dx,dy)),
            "depth_from_top_A": float(depth_A),
            "source_center_xy_A": source_xy.tolist(),
        },
        "source_neighborhood_R30_skin20": {
            "voxel_count": int(source_neighborhood.sum()),
            "Te_mean_K": float(te[source_neighborhood].mean()),
            "Te_max_K": float(te[source_neighborhood].max()),
            "Tl_mean_K": float(tl[source_neighborhood].mean()),
            "Tl_max_K": float(tl[source_neighborhood].max()),
            "definition": "Ag voxels within lateral r<=30 A of the source center (periodic x/y) and depth<=20 A from the top",
        },
        "high_temperature_regions": {
            "Te_50pct_peak_excess": zone_metrics(te,0.5),
            "Te_90pct_peak_excess": zone_metrics(te,0.9),
            "Tl_50pct_peak_excess": zone_metrics(tl,0.5),
            "Tl_90pct_peak_excess": zone_metrics(tl,0.9),
        },
        "phase": {"ag_liquid_volume_fraction": float(fl[ag].mean()), "max_liquid_fraction": float(fl[ag].max()),
                  "liquid_voxels": int(np.count_nonzero(fl[ag] > 0))},
        "energy_accounting": {
            "source_J": source_total,
            "electron_energy_J": float(final["electron_energy_J"]),
            "lattice_enthalpy_J": float(final["lattice_enthalpy_J"]),
            "electron_to_lattice_cumulative_J": float(final["electron_to_lattice_cumulative_J"]),
            "electron_diffusion_net_J_logged": float(final["electron_diffusion_net_J"]),
            "electron_diffusion_net_J_reconstructed": diffusion_net_J,
            "local_electron_redistribution_total_variation_J": diffusion_gross_J,
            "instantaneous_final_electron_face_flux_abs_W": face_flux_abs_W,
            "instantaneous_final_electron_face_flux_abs_by_axis_W_xyz": axis_gross_W,
            "instantaneous_final_net_power_out_of_te_hot_region_W": hot_net_out_W,
            "energy_balance_residual_J": float(final["energy_balance_residual_J"]),
            "energy_residual_relative": float(final["energy_residual_relative_to_source"]),
            "predicted_mean_dT_K_if_all_1ps_source_instantly_to_lattice": float(predicted_mean_dT_all_to_lattice),
            "predicted_mean_dT_K_from_current_lattice_enthalpy": float(predicted_mean_dT_coupled),
        },
        "electron_transport_at_Te_ag_mean": {
            "Te_mean_K": temp_at_mean, "Ce_J_m3_K": Ce, "ke_W_m_K": K, "G_W_m3_K": G,
            "alpha_m2_s": alpha, "sqrt_alpha_t_nm_at_0p1_0p5_1ps": diffusion_lengths_nm.tolist(),
            "Ce_over_G_ps": float(electron_tau_ps), "rhoCp_over_G_ps": float(lattice_tau_ps),
            "property_table_range_K": [float(temp_table.min()),float(temp_table.max())],
        },
        "time_samples": samples,
        "limitations": [
            "Intermediate scalar histories are available every 0.025 ps; full volumetric fields survive only at 0.15 ps and 1 ps in this checkpoint chain.",
            "The final instantaneous electron face flux is reconstructed from the final temperature field and binary Ag mask; transient flux history was not saved.",
            "Boundary A treats Ti3SiC2 as a geometric mask only; it has no heat capacity or thermal conductivity in this run. z is insulated.",
            "Electronic parameters are first-pass candidate data, 300–20000 K, with cross-source provenance requiring independent verification.",
        ],
    }
    (args.out_dir / "audit_metrics.json").write_text(json.dumps(metrics,indent=2)+"\n")

    # Plots: fields masked to true Ag; z slice one skin depth below top.
    zidx = int(np.clip(round(ag.shape[2]-1 - 20.0/spacing_A[2]),0,ag.shape[2]-1))
    yidx = ag.shape[1]//2
    fig, axes = plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
    panels=[(te[:,:,zidx].T,ag[:,:,zidx].T,"Electron temperature $T_e$ (K)","xy"),
            (tl[:,:,zidx].T,ag[:,:,zidx].T,"Lattice temperature $T_l$ (K)","xy"),
            (te[:,yidx,:].T,ag[:,yidx,:].T,"Electron temperature $T_e$ (K)","xz"),
            (tl[:,yidx,:].T,ag[:,yidx,:].T,"Lattice temperature $T_l$ (K)","xz")]
    for ax,(a,m,title,plane) in zip(axes.flat,panels):
        masked=np.ma.masked_where(~m,a)
        im=ax.imshow(masked,origin="lower",aspect="auto",interpolation="nearest",extent=(0,size_A[0]/10,0,(size_A[1] if plane=="xy" else size_A[2])/10))
        ax.set_title(title+f"; {plane} slice")
        ax.set_xlabel("x (nm)"); ax.set_ylabel("y (nm)" if plane=="xy" else "z (nm)")
        fig.colorbar(im,ax=ax,label="K")
    fig.suptitle(f"QSGS native Ag graph, t=1 ps; z(xy)={zidx}, y(xz)={yidx}")
    fig.savefig(args.out_dir/"temperature_fields_1ps_xy_xz.png",dpi=180)
    plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
    tt=np.array([r["time_ps"] for r in progress["samples"]])
    for ax, names, title in [(axes[0],["electron_temperature_max_K","electron_temperature_mean_K","lattice_temperature_max_K","lattice_temperature_mean_K"],"Ag-weighted temperature history"),
                              (axes[1],["electron_energy_J","lattice_enthalpy_J","source_energy_J","electron_to_lattice_cumulative_J"],"Energy ledger")]:
        for key in names:
            ys=np.array([r.get(key,np.nan) for r in progress["samples"]])
            ax.plot(tt,ys,label=key.replace("_K",""))
        ax.set_title(title); ax.set_xlabel("time (ps)"); ax.grid(True,alpha=.3); ax.legend(fontsize=8)
    axes[1].set_ylabel("J")
    fig.savefig(args.out_dir/"temperature_energy_history_1ps.png",dpi=180)
    plt.close(fig)

    # Reconstructed net local redistribution: positive means cumulative net
    # electron energy gained from diffusion; negative means lost by diffusion.
    sl=diffusion_J_cell[:,yidx,:].T/EV_J
    mm=ag[:,yidx,:].T
    lim=float(np.quantile(np.abs(sl[mm]),0.99)) if np.any(mm) else 1.0
    fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
    im=ax.imshow(np.ma.masked_where(~mm,sl),origin="lower",aspect="auto",cmap="coolwarm",vmin=-lim,vmax=lim,
                 extent=(0,size_A[0]/10,0,size_A[2]/10),interpolation="nearest")
    ax.set_xlabel("x (nm)"); ax.set_ylabel("z (nm)"); ax.set_title("Reconstructed cumulative local electron diffusion energy per Ag voxel (eV), y midplane")
    fig.colorbar(im,ax=ax,label="eV/voxel; red=net gain, blue=net loss")
    fig.savefig(args.out_dir/"electron_local_diffusion_redistribution_1ps.png",dpi=180)
    plt.close(fig)
    print(json.dumps({"metrics":str(args.out_dir/"audit_metrics.json"),"plots":["temperature_fields_1ps_xy_xz.png","temperature_energy_history_1ps.png","electron_local_diffusion_redistribution_1ps.png"],"te":te_stats,"tl":tl_stats,"lattice_max_location":metrics["lattice_max_location"],"energy":metrics["energy_accounting"],"transport":metrics["electron_transport_at_Te_ag_mean"]},indent=2))


if __name__ == "__main__":
    main()
