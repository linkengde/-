#!/usr/bin/env python3
"""Geometry-only stochastic bicontinuous Ag/Ti3SiC2 reconstruction.

This makes a reproducible spinodal-like 3D morphology approximation from a
filtered random field, then populates it with independently generated ideal
Ag and Ti3SiC2 parent lattices. It does not read atom coordinates from a melt
trajectory or reuse the layered heating model.
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
from ase.build import make_supercell
from ase.io import read, write
from ase.neighborlist import neighbor_list
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


ROOT = Path(__file__).resolve().parent
REPO = Path.cwd()
BASE = REPO / "research/mace-v12-transfer/coordination/reports/window-f/bicontinuous_ag_ti3sic2_ttm_20261015"
SOURCE = BASE / "exploratory_melt_cool_20261018/input/NEW_MODEL_OVITO.extxyz"
CIF = REPO / "research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif"
SEED = 20261021
TARGET_AG_WT = 0.5995
FIELD_SIGMA_GRID = (5.0, 9.0, 9.0)
GRID_SHAPE = (96, 160, 160)
AG_A0 = 4.09


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parent_sites():
    ref = read(SOURCE)
    cell_ref = np.asarray(ref.cell.array, float)
    x_lo, x_hi = float(ref.positions[:, 0].min()), float(ref.positions[:, 0].max())
    ly, lz = float(cell_ref[1, 1]), float(cell_ref[2, 2])
    # Reconstruct Ti3SiC2 from COD 9009647. c || X, a1 || Y, a1+2a2 || Z.
    tsc0 = read(CIF)
    if tsc0.get_chemical_formula() != "C4Si2Ti6":
        raise RuntimeError("Unexpected Ti3SiC2 CIF formula")
    pmat = np.array([[0, 0, 1], [1, 0, 0], [1, 2, 0]], int)
    tsc_cell = make_supercell(tsc0, pmat)
    rotate = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], float)  # (x,y,z)->(z,x,y)
    tsc_cell.positions = tsc_cell.positions @ rotate
    tsc_cell.cell.array = tsc_cell.cell.array @ rotate
    tsc_cell.wrap()
    ax, ay, az = np.asarray(tsc_cell.cell.lengths(), float)
    ny, nz = int(round(ly / ay)), int(round(lz / az))
    strain_y, strain_z = ly / (ny * ay), lz / (nz * az)
    if abs(strain_y - 1.0) > 0.01 or abs(strain_z - 1.0) > 0.01:
        raise RuntimeError("Ti3SiC2 lateral matching strain exceeds 1 percent")
    tsc_pos, tsc_sym = [], []
    for ix in range(int(np.ceil((x_hi - x_lo) / ax))):
        for iy in range(ny):
            for iz in range(nz):
                q = tsc_cell.positions.copy()
                q[:, 0] += x_lo + ix * ax
                q[:, 1] = (q[:, 1] + iy * ay) * strain_y
                q[:, 2] = (q[:, 2] + iz * az) * strain_z
                keep = (q[:, 0] <= x_hi) & (q[:, 1] < ly) & (q[:, 2] < lz)
                tsc_pos.append(q[keep])
                tsc_sym.extend(np.asarray(tsc_cell.get_chemical_symbols())[keep].tolist())
    tsc_pos = np.vstack(tsc_pos)

    # Orthorhombic primitive FCC cell: X=[001], Y=[110], Z=[1-10].
    ag_ay, ag_az = AG_A0 / np.sqrt(2.0), AG_A0 / np.sqrt(2.0)
    ag_basis = np.array([[0, 0, 0], [AG_A0 / 2.0, ag_ay / 2.0, ag_az / 2.0]])
    ax_ag = AG_A0
    nxx = int(np.ceil((x_hi - x_lo) / ax_ag))
    ny_ag, nz_ag = int(round(ly / ag_ay)), int(round(lz / ag_az))
    strain_ag_y, strain_ag_z = ly / (ny_ag * ag_ay), lz / (nz_ag * ag_az)
    if abs(strain_ag_y - 1.0) > 0.01 or abs(strain_ag_z - 1.0) > 0.01:
        raise RuntimeError("Ag lateral matching strain exceeds 1 percent")
    ix, iy, iz, ib = np.meshgrid(np.arange(nxx), np.arange(ny_ag), np.arange(nz_ag), np.arange(2), indexing="ij")
    idx = np.column_stack((ix.ravel(), iy.ravel(), iz.ravel(), ib.ravel()))
    ag_pos = ag_basis[idx[:, 3]].copy()
    ag_pos[:, 0] += x_lo + idx[:, 0] * ax_ag
    ag_pos[:, 1] = (ag_pos[:, 1] + idx[:, 1] * ag_ay) * strain_ag_y
    ag_pos[:, 2] = (ag_pos[:, 2] + idx[:, 2] * ag_az) * strain_ag_z
    keep = (ag_pos[:, 0] <= x_hi) & (ag_pos[:, 1] < ly) & (ag_pos[:, 2] < lz)
    ag_pos = ag_pos[keep]
    cell = np.diag([cell_ref[0, 0], ly, lz])
    return ref, x_lo, x_hi, ly, lz, cell, tsc_pos, np.asarray(tsc_sym), ag_pos, (strain_y, strain_z, strain_ag_y, strain_ag_z)


def generate_field():
    rng = np.random.default_rng(SEED)
    field = rng.standard_normal(GRID_SHAPE)
    field = gaussian_filter(field, sigma=FIELD_SIGMA_GRID, mode=("reflect", "wrap", "wrap"))
    field -= field.mean()
    field /= field.std()
    return field


def sample_field(field, xlo, xhi, ly, lz, points):
    coords = np.vstack((
        (points[:, 0] - xlo) / (xhi - xlo) * (GRID_SHAPE[0] - 1),
        (points[:, 1] % ly) / ly * GRID_SHAPE[1],
        (points[:, 2] % lz) / lz * GRID_SHAPE[2],
    ))
    sampled = map_coordinates(field, coords, order=1, mode="wrap", prefilter=False)
    return sampled


def choose_threshold(ag_scores, tsc_scores):
    m_ag = 107.8682
    m_ti, m_si, m_c = 47.867, 28.085, 12.011
    tsc_mass = {"Ti": m_ti, "Si": m_si, "C": m_c}
    tsc_symbols = tsc_scores[1]
    lo = min(float(ag_scores.min()), float(tsc_scores[0].min())) - 1e-6
    hi = max(float(ag_scores.max()), float(tsc_scores[0].max())) + 1e-6
    for _ in range(48):
        mid = (lo + hi) / 2
        nag = int(np.count_nonzero(ag_scores > mid))
        tkeep = tsc_scores[0] <= mid
        mt = sum(np.count_nonzero(tkeep & (tsc_symbols == el)) * mass for el, mass in tsc_mass.items())
        frac = nag * m_ag / (nag * m_ag + mt)
        if frac > TARGET_AG_WT:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def phase_graph(atoms, phases, label, cutoff):
    i, j = neighbor_list("ij", atoms, cutoff)
    keep = (phases[i] == label) & (phases[j] == label) & (i < j)
    nodes = np.flatnonzero(phases == label)
    lut = np.full(len(atoms), -1, int)
    lut[nodes] = np.arange(len(nodes))
    u, v = lut[i[keep]], lut[j[keep]]
    mat = coo_matrix((np.ones(2 * len(u), np.uint8), (np.r_[u, v], np.r_[v, u])), shape=(len(nodes), len(nodes))).tocsr()
    ncomp, labs = connected_components(mat, directed=False)
    sizes = np.bincount(labs) if len(labs) else np.array([], int)
    return {"phase": label, "cutoff_A": cutoff, "components": int(ncomp),
            "largest_component_atoms": int(sizes.max()) if len(sizes) else 0,
            "largest_component_fraction": float(sizes.max() / len(nodes)) if len(nodes) else 0.0}


def main():
    OUT = ROOT
    OUT.mkdir(parents=True, exist_ok=True)
    ref, xlo, xhi, ly, lz, cell, tsc_pos, tsc_sym, ag_pos, strains = parent_sites()
    grid = generate_field()
    ag_scores = sample_field(grid, xlo, xhi, ly, lz, ag_pos)
    tsc_scores = sample_field(grid, xlo, xhi, ly, lz, tsc_pos)
    threshold = choose_threshold(ag_scores, (tsc_scores, tsc_sym))
    ag_keep = ag_scores > threshold
    tsc_keep = tsc_scores <= threshold
    positions = np.vstack((tsc_pos[tsc_keep], ag_pos[ag_keep]))
    symbols = np.r_[tsc_sym[tsc_keep], np.full(int(ag_keep.sum()), "Ag")]
    phases = np.r_[np.full(int(tsc_keep.sum()), "Ti3SiC2"), np.full(int(ag_keep.sum()), "Ag")]
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=[False, True, True])
    ids = np.arange(1, len(atoms) + 1, dtype=int)
    atoms.arrays["atom_id"] = ids
    atoms.arrays["phase"] = phases.astype("U8")
    atoms.arrays["random_field_value"] = np.r_[tsc_scores[tsc_keep], ag_scores[ag_keep]]
    atoms.arrays["region"] = np.where(phases == "Ag", "Ag_network", "TSC_network").astype("U16")
    atoms.info.update({"model_status": "GEOMETRY_ONLY_RANDOM_BICONTINUOUS_APPROXIMATION_NOT_FORCE_VALIDATED",
                       "morphology": "filtered-Gaussian random-field spinodal-like approximation",
                       "random_seed": SEED, "field_sigma_grid": list(FIELD_SIGMA_GRID),
                       "field_threshold": float(threshold), "blue_display_phase": "Ag",
                       "source_reference_geometry": str(SOURCE.relative_to(REPO))})

    # Geometry and interface audit.
    pairs_i, pairs_j, pairs_d = neighbor_list("ijd", atoms, 4.5)
    cross = (phases[pairs_i] != phases[pairs_j]) & (pairs_i < pairs_j)
    cross_dist = pairs_d[cross]
    graphs = [phase_graph(atoms, phases, "Ag", 3.15), phase_graph(atoms, phases, "Ti3SiC2", 2.75)]
    counts = {el: int(np.count_nonzero(symbols == el)) for el in ("Ag", "Ti", "Si", "C")}
    m = {"Ag": 107.8682, "Ti": 47.867, "Si": 28.085, "C": 12.011}
    ag_mass = counts["Ag"] * m["Ag"]
    tsc_mass = sum(counts[el] * m[el] for el in ("Ti", "Si", "C"))
    audit = {
        "model_status": "GEOMETRY_ONLY_NOT_FORCE_VALIDATED",
        "requested_reference": "user_uploaded 2D irregular bicontinuous pattern; 3D is necessarily a stochastic approximation",
        "construction": "independent ideal FCC Ag + COD 9009647 Ti3SiC2 parent crystals classified by a filtered-Gaussian random field; no trajectory coordinates used",
        "seed": SEED, "field_sigma_grid": list(FIELD_SIGMA_GRID), "grid_shape": list(GRID_SHAPE), "field_threshold": float(threshold),
        "atoms": int(len(atoms)), "element_counts": counts,
        "Ag_mass_fraction": float(ag_mass / (ag_mass + tsc_mass)),
        "box_A": cell.diagonal().tolist(), "material_x_extent_A": [xlo, xhi], "pbc": [False, True, True],
        "lateral_strains": {"TSC_Y": strains[0] - 1, "TSC_Z": strains[1] - 1, "Ag_Y": strains[2] - 1, "Ag_Z": strains[3] - 1},
        "phase_graphs": graphs,
        "cross_phase_min_A": float(cross_dist.min()) if len(cross_dist) else None,
        "cross_phase_percentiles_A": np.percentile(cross_dist, [1, 5, 50, 95]).tolist() if len(cross_dist) else [],
        "cross_phase_pairs_lt_2_A": int(np.count_nonzero(cross_dist < 2.0)),
        "source_cif_sha256": sha(CIF), "reference_extxyz_sha256": sha(SOURCE),
        "crystal_generation_strains_note": "ideal lattices are independently built; only small common Y/Z matching strain is applied",
        "scientific_limit": "single uploaded 2D image does not determine a unique 3D morphology; this is one reproducible spinodal-like analog, not an exact reconstruction",
    }

    write(OUT / "RANDOM_BICONTINUOUS.extxyz", atoms, format="extxyz")
    write(OUT / "RANDOM_BICONTINUOUS.data", atoms, format="lammps-data", atom_style="atomic", specorder=["Ti", "Si", "C", "Ag"])
    with (OUT / "PHASE_ATOM_GROUPS.csv").open("w") as f:
        f.write("atom_id,element,phase,region,field_value\n")
        for i, (sym, phase, region, value) in enumerate(zip(symbols, phases, atoms.arrays["region"], atoms.arrays["random_field_value"]), 1):
            f.write(f"{i},{sym},{phase},{region},{value:.8f}\n")
    (OUT / "GEOMETRY_AUDIT.json").write_text(json.dumps(audit, indent=2) + "\n")
    render_slices(grid, threshold, xlo, xhi, ly, lz, OUT / "RANDOM_BICONTINUOUS_SLICES.png")
    return audit


def render_slices(field, threshold, xlo, xhi, ly, lz, path):
    ny, nz = field.shape[1:]
    cmap = plt.matplotlib.colors.ListedColormap(["#d5d7dc", "#2864d2"])
    xs = [0.15, 0.35, 0.55, 0.75, 0.95]
    fig, axes = plt.subplots(1, 5, figsize=(14, 3.4), dpi=170)
    for ax, frac in zip(axes, xs):
        ix = min(int(frac * (field.shape[0] - 1)), field.shape[0] - 1)
        phase = field[ix] > threshold
        ax.imshow(phase.T, origin="lower", extent=(0, ly, 0, lz), interpolation="nearest",
                  cmap=cmap, vmin=0, vmax=1, aspect="equal")
        ax.set_title(f"X = {xlo + frac * (xhi - xlo):.1f} Å")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("Spinodal-like 3D random-field approximation: YZ slices\nBlue = Ag (provisional); gray = Ti₃SiC₂", weight="bold")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
