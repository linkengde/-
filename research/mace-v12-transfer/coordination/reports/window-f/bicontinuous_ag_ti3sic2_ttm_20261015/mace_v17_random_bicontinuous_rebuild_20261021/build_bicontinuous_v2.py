#!/usr/bin/env python3
"""Build an independent, irregular 3-D Ag/Ti3SiC2 crystalline morphology.

This is a geometry reconstruction from ideal parent crystals and a seeded
filtered random field. It does not read an MD trajectory. Interface overlap
resolution is deterministic and recorded; the result is geometry-only and is
not a force-field validation or a heating-ready declaration.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.io import write
from ase.neighborlist import neighbor_list
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

import build_random_bicontinuous as base


ROOT = Path(__file__).resolve().parent
SEED = 20261021
base.SEED = SEED
base.FIELD_SIGMA_GRID = (6.0, 9.0, 9.0)
base.GRID_SHAPE = (128, 160, 160)
base.TARGET_AG_WT = 0.608  # compensate composition change from explicit interface-site exclusion
TARGET_AG_WT = 0.5995
MIN_CROSS_DISTANCE_A = 2.20


def periodic_graph(atoms: Atoms, phases: np.ndarray, phase: str, cutoff: float) -> dict:
    """Return component size and periodic winding information for one phase."""
    i, j, shifts = neighbor_list("ijS", atoms, cutoff)
    members = np.flatnonzero(phases == phase)
    local = np.full(len(atoms), -1, dtype=int)
    local[members] = np.arange(len(members))
    keep = (phases[i] == phase) & (phases[j] == phase)
    ii, jj, ss = local[i[keep]], local[j[keep]], shifts[keep].astype(int)
    n = len(members)
    rows = np.r_[ii, jj]
    cols = np.r_[jj, ii]
    mat = coo_matrix((np.ones(len(rows), np.uint8), (rows, cols)), shape=(n, n)).tocsr()
    ncomp, labels = connected_components(mat, directed=False)
    sizes = np.bincount(labels) if n else np.array([], dtype=int)

    # Lift each component through the periodic graph. A non-zero cycle residual
    # proves the phase winds through that periodic cell axis.
    adjacency: list[list[tuple[int, np.ndarray]]] = [[] for _ in range(n)]
    for u, v, s in zip(ii, jj, ss):
        adjacency[int(u)].append((int(v), s))
        adjacency[int(v)].append((int(u), -s))
    seen = np.zeros(n, dtype=bool)
    image = np.zeros((n, 3), dtype=int)
    winding = np.zeros(3, dtype=bool)
    for root in range(n):
        if seen[root]:
            continue
        seen[root] = True
        stack = [root]
        while stack:
            u = stack.pop()
            for v, shift in adjacency[u]:
                proposed = image[u] + shift
                if not seen[v]:
                    seen[v] = True
                    image[v] = proposed
                    stack.append(v)
                else:
                    winding |= (proposed != image[v])
    largest = int(np.argmax(sizes)) if n else -1
    largest_nodes = members[labels == largest] if n else np.array([], dtype=int)
    return {
        "phase": phase,
        "neighbor_cutoff_A": cutoff,
        "components": int(ncomp),
        "largest_component_atoms": int(sizes.max()) if n else 0,
        "largest_component_fraction": float(sizes.max() / n) if n else 0.0,
        "largest_component_X_extent_A": [float(atoms.positions[largest_nodes, 0].min()),
                                         float(atoms.positions[largest_nodes, 0].max())] if len(largest_nodes) else [],
        "periodic_winding": {"X_nonperiodic": False, "Y": bool(winding[1]), "Z": bool(winding[2])},
    }


def remove_only_overlapping_interface_sites(atoms: Atoms, phases: np.ndarray,
                                            field_values: np.ndarray, threshold: float) -> tuple[Atoms, np.ndarray, dict]:
    """Resolve unlike-phase clashes with a reproducible boundary-site rule.

    For each unlike pair closer than MIN_CROSS_DISTANCE_A, retain the site
    deeper on its assigned side of the phase field and remove the boundary-side
    site. This is an interface construction operation, not post hoc force
    cleanup. Every deleted ID and its original phase/field score is reported.
    """
    alive = np.ones(len(atoms), dtype=bool)
    depth = np.abs(field_values - threshold)
    removed: list[int] = []
    initial_short = 0
    while True:
        current = np.flatnonzero(alive)
        a = atoms[current]
        p = phases[current]
        ii, jj, dd = neighbor_list("ijd", a, MIN_CROSS_DISTANCE_A)
        cross = (p[ii] != p[jj]) & (ii < jj)
        ii, jj, dd = ii[cross], jj[cross], dd[cross]
        if len(dd) == 0:
            break
        if not initial_short:
            initial_short = int(len(dd))
        # Resolve the tightest contacts first; ties remove the site with more
        # conflict partners, then deterministically prefer the Ag site.
        order = np.argsort(dd, kind="stable")
        degree = np.bincount(np.r_[ii, jj], minlength=len(a))
        removed_this_pass = 0
        for edge in order:
            u, v = int(ii[edge]), int(jj[edge])
            gu, gv = int(current[u]), int(current[v])
            if not alive[gu] or not alive[gv]:
                continue
            du, dv = depth[gu], depth[gv]
            if du < dv - 1e-10:
                victim = gu
            elif dv < du - 1e-10:
                victim = gv
            elif degree[u] != degree[v]:
                victim = gu if degree[u] > degree[v] else gv
            else:
                # Fixed tie break makes the file byte-for-byte reproducible.
                victim = gu if phases[gu] == "Ag" else gv
            alive[victim] = False
            removed.append(victim)
            removed_this_pass += 1
        if not removed_this_pass:
            raise RuntimeError("Could not resolve cross-phase overlap graph")
    out = atoms[alive]
    out_phases = phases[alive]
    return out, out_phases, {
        "minimum_unlike_phase_distance_A": MIN_CROSS_DISTANCE_A,
        "initial_unlike_pairs_below_minimum": initial_short,
        "removed_interface_sites": len(removed),
        "removed_by_phase": {
            phase: int(sum(phases[idx] == phase for idx in removed))
            for phase in ("Ag", "Ti3SiC2")
        },
        "removed_original_zero_based_indices": removed,
        "survivors": int(np.count_nonzero(alive)),
    }


def render_views(atoms: Atoms, phases: np.ndarray, field: np.ndarray,
                 threshold: float, xlo: float, xhi: float, ly: float, lz: float) -> None:
    # Orthographic projection looking from the positive X/Y/Z corner.
    pos = atoms.positions
    view = np.array([2.3, -1.1, 1.0])
    view /= np.linalg.norm(view)
    up0 = np.array([0.0, 0.0, 1.0])
    right = np.cross(up0, view)
    right /= np.linalg.norm(right)
    up = np.cross(view, right)
    xy = np.column_stack((pos @ right, pos @ up))
    depth = pos @ view
    order = np.argsort(depth)
    palette = {"Ag": "#f0a320", "Ti3SiC2": "#d83dc3"}
    colors = np.array([palette[p] for p in phases])
    fig, ax = plt.subplots(figsize=(8.8, 6.7), dpi=190)
    ax.scatter(xy[order, 0], xy[order, 1], c=colors[order], s=1.35,
               linewidths=0, rasterized=True, alpha=0.94)
    corners = np.array([[x, y, z] for x in (xlo, xhi) for y in (0.0, ly) for z in (0.0, lz)])
    corner_xy = np.column_stack((corners @ right, corners @ up))
    for a in range(8):
        for b in range(a + 1, 8):
            if sum(corners[a, k] != corners[b, k] for k in range(3)) == 1:
                ax.plot(corner_xy[[a, b], 0], corner_xy[[a, b], 1], color="#333333", lw=0.45, alpha=0.55)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("3-D bicontinuous unit cell · Ag (gold) / Ti₃SiC₂ (magenta)\n"
                 "Oblique view; Y/Z edges are periodic", fontsize=11)
    fig.tight_layout(pad=0.2)
    fig.savefig(ROOT / "BICONNECTIVE_3D_ATOM_VIEW.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), dpi=180)
    field_x = np.linspace(xlo, xhi, field.shape[0])
    xs = [0.18, 0.52, 0.86]
    cmap = matplotlib.colors.ListedColormap(["#d83dc3", "#f0a320"])
    for ax, frac in zip(axes, xs):
        ix = int(frac * (field.shape[0] - 1))
        # The pixel morphology is paired with the actual atom-level model view.
        ax.imshow((field[ix] > threshold).T, origin="lower", extent=(0, ly, 0, lz),
                  cmap=cmap, interpolation="nearest", aspect="equal", vmin=0, vmax=1)
        ax.set_title(f"YZ slice, X={field_x[ix]:.1f} Å")
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Three distinct sections through the same random 3-D phase field\n"
                 "gold = Ag · magenta = Ti₃SiC₂", fontsize=11)
    fig.tight_layout()
    fig.savefig(ROOT / "BICONNECTIVE_YZ_SLICES.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    ref, xlo, xhi, ly, lz, cell, tsc_pos, tsc_sym, ag_pos, strains = base.parent_sites()
    field = base.generate_field()

    # A broad, smooth Ag enrichment near the central free-surface region gives
    # a candidate local heating site without turning the entire top into a slab.
    gx = np.linspace(0.0, 1.0, field.shape[0])[:, None, None]
    gy = np.arange(field.shape[1])[None, :, None] / field.shape[1]
    gz = np.arange(field.shape[2])[None, None, :] / field.shape[2]
    # Wrapped lateral distances; Gaussian width is about 14 Å in Y/Z.
    dy = np.minimum(np.abs(gy - 0.5), 1.0 - np.abs(gy - 0.5)) * ly
    dz = np.minimum(np.abs(gz - 0.5), 1.0 - np.abs(gz - 0.5)) * lz
    xdist = (1.0 - gx) * (xhi - xlo)
    field += 1.1 * np.exp(-0.5 * ((xdist / 10.0) ** 2 + (dy / 14.0) ** 2 + (dz / 14.0) ** 2))
    field -= field.mean()
    field /= field.std()

    ag_score = base.sample_field(field, xlo, xhi, ly, lz, ag_pos)
    tsc_score = base.sample_field(field, xlo, xhi, ly, lz, tsc_pos)
    threshold = base.choose_threshold(ag_score, (tsc_score, tsc_sym))
    ag_keep = ag_score > threshold
    tsc_keep = tsc_score <= threshold
    pos = np.vstack((tsc_pos[tsc_keep], ag_pos[ag_keep]))
    syms = np.r_[tsc_sym[tsc_keep], np.full(int(ag_keep.sum()), "Ag")]
    phases = np.r_[np.full(int(tsc_keep.sum()), "Ti3SiC2"), np.full(int(ag_keep.sum()), "Ag")]
    scores = np.r_[tsc_score[tsc_keep], ag_score[ag_keep]]
    atoms = Atoms(symbols=syms, positions=pos, cell=cell, pbc=[False, True, True])
    assigned_phases, assigned_symbols = phases.copy(), syms.copy()
    assigned_scores, assigned_positions = scores.copy(), pos.copy()
    atoms, phases, overlap_audit = remove_only_overlapping_interface_sites(atoms, phases, scores, threshold)
    kept_original = np.ones(len(assigned_phases), dtype=bool)
    removed_indices = overlap_audit.pop("removed_original_zero_based_indices")
    kept_original[removed_indices] = False
    scores = assigned_scores[kept_original]
    with (ROOT / "INTERFACE_SITE_REMOVALS.csv").open("w") as f:
        f.write("source_index_0based,element,phase,x_A,y_A,z_A,field_value\n")
        for k in removed_indices:
            x, y, z = assigned_positions[k]
            f.write(f"{k},{assigned_symbols[k]},{assigned_phases[k]},{x:.8f},{y:.8f},{z:.8f},{assigned_scores[k]:.8f}\n")

    ids = np.arange(1, len(atoms) + 1, dtype=int)
    atoms.arrays["atom_id"] = ids
    atoms.arrays["phase"] = phases.astype("U8")
    atoms.arrays["field_value"] = scores
    atoms.arrays["interface_site"] = (np.abs(scores - threshold) < 0.20)
    atoms.arrays["region"] = np.where(
        (phases == "Ag") & (atoms.positions[:, 0] > xhi - 15.0)
        & (np.linalg.norm((atoms.positions[:, 1:] - np.array([ly / 2, lz / 2]) + np.array([ly / 2, lz / 2])) % np.array([ly, lz]) - np.array([ly / 2, lz / 2]), axis=1) < 18.0),
        "Ag_surface_enriched", np.where(phases == "Ag", "Ag_network", "TSC_backbone"),
    ).astype("U24")

    # Cross-phase distances after boundary construction.
    ii, jj, dd = neighbor_list("ijd", atoms, 4.5)
    cross = (phases[ii] != phases[jj]) & (ii < jj)
    cross_d = dd[cross]
    graphs = [periodic_graph(atoms, phases, "Ag", 3.15),
              periodic_graph(atoms, phases, "Ti3SiC2", 2.75)]
    counts = {el: int(np.count_nonzero(np.array(atoms.get_chemical_symbols()) == el)) for el in ("Ag", "Ti", "Si", "C")}
    masses = {"Ag": 107.8682, "Ti": 47.867, "Si": 28.085, "C": 12.011}
    ag_mass = counts["Ag"] * masses["Ag"]
    tsc_mass = sum(counts[e] * masses[e] for e in ("Ti", "Si", "C"))
    domain_vol = float(np.mean(field > threshold))
    phase_vol = {"Ag": domain_vol, "Ti3SiC2": 1.0 - domain_vol}
    graph_sensitivity = {
        "Ag": [periodic_graph(atoms, phases, "Ag", c) for c in (2.95, 3.15, 3.30)],
        "Ti3SiC2": [periodic_graph(atoms, phases, "Ti3SiC2", c) for c in (2.55, 2.75, 2.90)],
    }
    audit = {
        "model_status": "GEOMETRY_ONLY_NOT_FORCE_VALIDATED_NOT_READY_FOR_HEATING",
        "reference_interpretation": "irregular 3-D atomistic block with interspersed domains; 2-D image does not define a unique 3-D volume",
        "construction": "independent ideal FCC Ag and COD 9009647 Ti3SiC2 lattices classified by a seeded filtered Gaussian field; no historical melt-pool coordinates used",
        "seed": SEED,
        "field_sigma_grid": list(base.FIELD_SIGMA_GRID), "field_grid_shape": list(base.GRID_SHAPE),
        "pre_interface_exclusion_target_Ag_mass_fraction": base.TARGET_AG_WT,
        "local_surface_Ag_bias": {"center": [xhi - 2.0, ly / 2, lz / 2], "sigma_A": [10.0, 14.0, 14.0], "amplitude": 1.1},
        "threshold": float(threshold), "number_of_atoms": int(len(atoms)), "element_counts": counts,
        "Ag_mass_fraction": float(ag_mass / (ag_mass + tsc_mass)),
        "phase_domain_volume_fraction_before_interface_site_exclusion": phase_vol,
        "cell_A": cell.diagonal().tolist(), "material_X_extent_A": [xlo, xhi], "pbc": [False, True, True],
        "lateral_strains": {"Ti3SiC2_Y": strains[0] - 1, "Ti3SiC2_Z": strains[1] - 1,
                            "Ag_Y": strains[2] - 1, "Ag_Z": strains[3] - 1},
        "interface_overlap_resolution": overlap_audit,
        "unlike_phase_distance_A": {"minimum": float(cross_d.min()) if len(cross_d) else None,
                                     "percentiles_1_5_50_95": np.percentile(cross_d, [1, 5, 50, 95]).tolist() if len(cross_d) else [],
                                     "pairs_below_2_A": int(np.count_nonzero(cross_d < 2.0)),
                                     "pairs_below_2_5_A": int(np.count_nonzero(cross_d < 2.5))},
        "phase_connectivity": graphs,
        "phase_connectivity_cutoff_sensitivity": graph_sensitivity,
        "Ag_surface_enriched_atom_count": int(np.count_nonzero(atoms.arrays["region"] == "Ag_surface_enriched")),
        "parent_CIF_sha256": base.sha(base.CIF), "reference_structure_sha256": base.sha(base.SOURCE),
        "interpretation_limit": "topology and geometry only; no force-field, relaxation, melting, or solidification validation",
    }
    (ROOT / "RANDOM_BICONTINUOUS_V2_AUDIT.json").write_text(json.dumps(audit, indent=2) + "\n")
    write(ROOT / "RANDOM_BICONTINUOUS_V2.extxyz", atoms, format="extxyz")
    write(ROOT / "RANDOM_BICONTINUOUS_V2.data", atoms, format="lammps-data", atom_style="atomic",
          specorder=["Ti", "Si", "C", "Ag"])
    with (ROOT / "PHASE_ATOM_GROUPS_V2.csv").open("w") as f:
        f.write("atom_id,element,phase,region,interface_site,field_value\n")
        for i, (sym, phase, region, is_iface, val) in enumerate(zip(
                atoms.get_chemical_symbols(), phases, atoms.arrays["region"],
                atoms.arrays["interface_site"], scores), 1):
            f.write(f"{i},{sym},{phase},{region},{int(is_iface)},{val:.8f}\n")
    render_views(atoms, phases, field, threshold, xlo, xhi, ly, lz)
    (ROOT / "THREE_D_BICONNECTIVITY_GEOMETRY_README.md").write_text(
        "# Irregular three-dimensional Ag–Ti₃SiC₂ bicontinuous geometry\n\n"
        "This is a reproducible approximation of the user's 3-D reference: both phases form interpenetrating, continuous volumes through the model. It uses independently generated ideal FCC Ag and ideal Ti₃SiC₂ parent crystals classified with a seeded smooth random field; it does not reuse coordinates from a melt trajectory and does not implement the referenced four-parameter generator. A single image cannot uniquely specify a 3-D volume.\n\n"
        "The output is a geometry candidate only. It has not been evaluated with MACE V17, relaxed, heated, or validated for force-field stability. The V1 candidate with sub-ångström overlaps is retained as a failed attempt and is not used.\n\n"
        f"- Atom count: {len(atoms):,}; Ag mass fraction: {ag_mass / (ag_mass + tsc_mass):.4%}.\n"
        f"- Box: {cell.diagonal()[0]:.3f} × {cell.diagonal()[1]:.3f} × {cell.diagonal()[2]:.3f} Å; PBC: X free, Y/Z periodic.\n"
        f"- Cross-phase nearest distance after deterministic interface-site construction: {cross_d.min():.3f} Å. This is a geometry check only, not proof of stable forces.\n"
        f"- Ag largest-component fraction: {graphs[0]['largest_component_fraction']:.5%}; X span: {graphs[0]['largest_component_X_extent_A'][0]:.2f}–{graphs[0]['largest_component_X_extent_A'][1]:.2f} Å; Y/Z winding: {graphs[0]['periodic_winding']['Y']}/{graphs[0]['periodic_winding']['Z']}.\n"
        f"- Ti₃SiC₂ largest-component fraction: {graphs[1]['largest_component_fraction']:.5%}; X span: {graphs[1]['largest_component_X_extent_A'][0]:.2f}–{graphs[1]['largest_component_X_extent_A'][1]:.2f} Å; Y/Z winding: {graphs[1]['periodic_winding']['Y']}/{graphs[1]['periodic_winding']['Z']}.\n"
        f"- Interface site selection removed {overlap_audit['removed_interface_sites']:,} boundary lattice sites to resolve initial unlike-phase contacts below {MIN_CROSS_DISTANCE_A:.2f} Å. Exact removed sites are listed in `INTERFACE_SITE_REMOVALS.csv`; this changes local interface coordination and requires subsequent static-force review before MD.\n\n"
        "Files: `RANDOM_BICONTINUOUS_V2.data` (LAMMPS atomic style), `RANDOM_BICONTINUOUS_V2.extxyz` (OVITO/ASE), `PHASE_ATOM_GROUPS_V2.csv`, `BICONNECTIVE_3D_ATOM_VIEW.png`, `BICONNECTIVE_YZ_SLICES.png`, and `RANDOM_BICONTINUOUS_V2_AUDIT.json`.\n"
    )
    print(json.dumps({k: audit[k] for k in ("number_of_atoms", "element_counts", "Ag_mass_fraction", "phase_connectivity", "unlike_phase_distance_A")}, indent=2))
    return audit


if __name__ == "__main__":
    main()
