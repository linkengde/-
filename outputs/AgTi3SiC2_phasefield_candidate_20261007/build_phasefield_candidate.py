#!/usr/bin/env python3
"""Build a phase-field Ag/Ti3SiC2 morphology and a coarse voxel proxy.

The phase-field parameters are morphology controls, not material-calibrated
thermodynamic parameters. No minimization or molecular dynamics is performed.
"""
from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
SOURCE = Path(__file__).resolve().parent / "inputs" / "planar_reference.data"
EXPECTED_SOURCE_SHA256 = "80b53f26f0665607753ab7a7677d1d504e99b0d75ef78ec70f6d6c82cca7615a"
GRID = (36, 36, 48)
SEEDS = range(440101, 440105)
CHECKPOINTS = (100, 125, 150, 180, 220, 280, 350, 420, 500)
DT = 0.002
MOBILITY = 10.0
KAPPA_A2 = 0.45
MEAN_ORDER_PARAMETER = 0.20  # positive phase TSC target ~= 60 vol.%
TARGET_TSC = 0.60
HEIGHT = 50.0
VACUUM = 20.0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_data(path: Path):
    bounds = {}
    tilt = (0.0, 0.0, 0.0)
    atoms = []
    in_atoms = False
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[-2:] == ["xlo", "xhi"]:
            bounds["x"] = tuple(map(float, fields[:2]))
        elif len(fields) >= 4 and fields[-2:] == ["ylo", "yhi"]:
            bounds["y"] = tuple(map(float, fields[:2]))
        elif len(fields) >= 4 and fields[-2:] == ["zlo", "zhi"]:
            bounds["z"] = tuple(map(float, fields[:2]))
        elif len(fields) >= 6 and fields[-3:] == ["xy", "xz", "yz"]:
            tilt = tuple(map(float, fields[:3]))
        if line.strip().startswith("Atoms"):
            in_atoms = True
            continue
        if in_atoms and len(fields) >= 5 and fields[0].isdigit() and fields[1].isdigit():
            atoms.append((int(fields[0]), int(fields[1]), *map(float, fields[2:5])))
    arr = np.asarray(atoms, dtype=float)
    if len(arr) != 4783:
        raise RuntimeError(f"Expected 4783 source atoms, found {len(arr)}")
    return bounds, tilt, arr


def frac_xy(x, y, lx, ly, xy):
    v = y / ly
    u = (x - xy * v) / lx
    return u % 1.0, v % 1.0


def topology(mask):
    """6-neighbor components, periodic in x/y and open in z."""
    nx, ny, nz = mask.shape
    tiled = np.tile(mask, (3, 3, 1))
    labels, _ = ndimage.label(tiled, ndimage.generate_binary_structure(3, 1))
    central = labels[nx:2 * nx, ny:2 * ny, :]
    vals, counts = np.unique(central[central > 0], return_counts=True)
    total = int(mask.sum())
    if not total:
        return {"voxels": 0, "components": 0, "largest_fraction": 0.0,
                "wrap_x": False, "wrap_y": False, "spans_z": False}
    lab = int(vals[int(np.argmax(counts))])
    base = labels[nx:2 * nx, ny:2 * ny, :]
    loc = np.argwhere(central == lab)
    return {
        "voxels": total,
        "components": int(len(vals)),
        "largest_fraction": float(counts.max() / total),
        "wrap_x": bool(np.any((base == lab) & (labels[2 * nx:3 * nx, ny:2 * ny, :] == lab))),
        "wrap_y": bool(np.any((base == lab) & (labels[nx:2 * nx, 2 * ny:3 * ny, :] == lab))),
        "spans_z": bool(loc[:, 2].min() == 0 and loc[:, 2].max() == nz - 1),
    }


def phasefield_spectrum(shape, lx, ly, xy, height):
    nx, ny, nz = shape
    # Fractional in-plane reciprocal modes in the triclinic 120-degree cell.
    ku = 2.0 * np.pi * np.fft.fftfreq(nx, d=1.0 / nx)[:, None, None]
    kv = 2.0 * np.pi * np.fft.fftfreq(ny, d=1.0 / ny)[None, :, None]
    # Inverse of the in-plane direct metric for a=(lx,0), b=(xy,ly).
    guu = 1.0 / (lx * lx) + xy * xy / (lx * lx * ly * ly)
    guv = -xy / (lx * ly * ly)
    gvv = 1.0 / (ly * ly)
    kxy2 = guu * ku * ku + 2.0 * guv * ku * kv + gvv * kv * kv
    # Even reflection in z imposes a no-flux condition at both free surfaces.
    kz = 2.0 * np.pi * np.fft.fftfreq(2 * nz, d=height / nz)[None, None, :]
    return kxy2 + kz * kz


def evolve(seed, shape, k2, stop_step):
    rng = np.random.default_rng(seed)
    raw = rng.standard_normal(shape)
    reflected_noise = np.concatenate((raw, raw[:, :, ::-1]), axis=2)
    spectral_noise = np.fft.fftn(reflected_noise)
    # Seed the linearly unstable spinodal band; all morphology is then evolved by Cahn-Hilliard.
    k2_peak = 1.0 / (2.0 * KAPPA_A2)
    band = (k2 >= 0.5 * k2_peak) & (k2 <= 1.5 * k2_peak)
    noise = np.fft.ifftn(spectral_noise * band).real[:, :, :shape[2]]
    noise -= noise.mean()
    noise /= max(noise.std(), 1e-12)
    phi = MEAN_ORDER_PARAMETER + 0.25 * noise
    phi -= phi.mean() - MEAN_ORDER_PARAMETER
    nx, ny, nz = shape
    k2_3d = k2
    denominator = 1.0 + MOBILITY * DT * KAPPA_A2 * k2_3d * k2_3d
    checkpoints = set(CHECKPOINTS)
    best = None
    field_snapshots = []
    for step in range(1, stop_step + 1):
        reflected = np.concatenate((phi, phi[:, :, ::-1]), axis=2)
        ph = np.fft.fftn(reflected)
        nonlinear = np.fft.fftn(reflected ** 3 - reflected)
        ph = (ph - MOBILITY * DT * k2_3d * nonlinear) / denominator
        phi = np.fft.ifftn(ph).real[:, :, :nz]
        if step in checkpoints:
            # Quantile threshold fixes the target phase fraction exactly.
            flat = phi.ravel()
            n_tsc = int(round(TARGET_TSC * flat.size))
            indices = np.argpartition(flat, flat.size - n_tsc)[flat.size - n_tsc:]
            mask = np.zeros(flat.size, dtype=bool)
            mask[indices] = True
            mask = mask.reshape(shape)
            tsc = topology(mask)
            ag = topology(~mask)
            topo_ok = float(np.std(phi)) >= 0.58 and all(
                d["largest_fraction"] >= 0.97 and d["wrap_x"] and d["wrap_y"] and d["spans_z"]
                for d in (tsc, ag)
            )
            # Prefer continuous phases with a larger smallest-phase backbone.
            score = min(tsc["largest_fraction"], ag["largest_fraction"]) + 0.01 * min(float(np.std(phi)), 1.0)
            candidate = (topo_ok, score, seed, step, phi.copy(), mask.copy(), tsc, ag)
            field_snapshots.append(candidate)
            if best is None or (candidate[0], candidate[1]) > (best[0], best[1]):
                best = candidate
            if topo_ok:
                return candidate, field_snapshots
    return best, field_snapshots


def crystal_lattices(source_atoms, lx, ly, xy):
    def layer_groups(types):
        layers = defaultdict(list)
        for row in source_atoms:
            if int(row[1]) in types:
                layers[round(float(row[4]), 6)].append(row)
        return dict(sorted(layers.items()))

    tsc_layers = layer_groups({1, 2, 3})
    signatures = {}
    for z, rows in tsc_layers.items():
        signatures[z] = tuple(sorted((int(r[1]), *[round(x, 6) for x in frac_xy(r[2], r[3], lx, ly, xy)]) for r in rows))
    z0 = min(tsc_layers)
    repeats = [z - z0 for z, sig in signatures.items() if z - z0 > 15.0 and sig == signatures[z0]]
    if not repeats:
        raise RuntimeError("Could not derive Ti3SiC2 c repeat from source")
    c_tsc = min(repeats)
    tsc_motif = [r.copy() for r in source_atoms if int(r[1]) in (1, 2, 3) and float(r[4]) < c_tsc - 1e-5]
    if len(tsc_motif) != 1728:
        raise RuntimeError(f"Unexpected TSC motif atom count: {len(tsc_motif)}")
    ag_layers = layer_groups({4})
    ag_z = list(ag_layers)
    d111 = float(np.median(np.diff(ag_z)))
    ag_motifs = [[r.copy() for r in ag_layers[z]] for z in ag_z[:3]]
    if any(len(layer) != 169 for layer in ag_motifs):
        raise RuntimeError("Unexpected Ag(111) layer population")
    ag_z0 = ag_z[0] - 2.0 * c_tsc
    if not 0.0 <= ag_z0 < d111:
        raise RuntimeError(f"Unexpected Ag layer phase {ag_z0}")
    return c_tsc, d111, ag_z0, tsc_motif, ag_motifs


def write_lammps_data(path, atoms, lx, ly, xy):
    masses = {1: (47.867, "Ti"), 2: (28.0855, "Si"), 3: (12.0107, "C"), 4: (107.8682, "Ag")}
    with path.open("w") as f:
        f.write("Phase-field Ag-Ti3SiC2 bicontinuous geometry; geometry-only, no relaxation\n\n")
        f.write(f"{len(atoms)} atoms\n4 atom types\n\n")
        f.write(f"0.0 {lx:.12f} xlo xhi\n0.0 {ly:.12f} ylo yhi\n0.0 {HEIGHT + VACUUM:.12f} zlo zhi\n")
        f.write(f"{xy:.12f} 0.0 0.0 xy xz yz\n\nMasses\n\n")
        for t, (mass, name) in masses.items():
            f.write(f"{t} {mass:.7f} # {name}\n")
        f.write("\nAtoms # atomic\n\n")
        for i, atom in enumerate(atoms, 1):
            f.write(f"{i} {atom['type']} {atom['x']:.8f} {atom['y']:.8f} {atom['z']:.8f}\n")


def write_voxel_proxy(path, mask, lx, ly, xy):
    """Write one coarse-grained point at each phase-field cell center for OVITO."""
    nx, ny, nz = mask.shape
    total = int(mask.size)
    with path.open("w") as f:
        f.write("Coarse phase-field voxel centers; visualization only, not atomistic/MD input\n\n")
        f.write(f"{total} atoms\n2 atom types\n\n")
        f.write(f"0.0 {lx:.12f} xlo xhi\n0.0 {ly:.12f} ylo yhi\n0.0 {HEIGHT + VACUUM:.12f} zlo zhi\n")
        f.write(f"{xy:.12f} 0.0 0.0 xy xz yz\n\nMasses\n\n1 1.0 # TSC_voxel\n2 1.0 # Ag_voxel\n\nAtoms # atomic\n\n")
        atom_id = 0
        for i, j, k in np.ndindex(mask.shape):
            atom_id += 1
            u, v = (i + 0.5) / nx, (j + 0.5) / ny
            x, y, z = u * lx + v * xy, v * ly, (k + 0.5) * HEIGHT / nz
            typ = 1 if mask[i, j, k] else 2
            f.write(f"{atom_id} {typ} {x:.8f} {y:.8f} {z:.8f}\n")


def min_cross_distances(tsc_xyz, ag_xyz, lx, ly, xy):
    a = np.array([lx, 0.0, 0.0])
    b = np.array([xy, ly, 0.0])
    images = np.concatenate([ag_xyz + i * a + j * b for i in (-1, 0, 1) for j in (-1, 0, 1)], axis=0)
    d, idx = cKDTree(images).query(tsc_xyz, k=1)
    return d, idx


def main():
    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)
    source_hash = sha256(SOURCE)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise RuntimeError(f"Source checksum mismatch: {source_hash}")
    bounds, tilt, source_atoms = parse_data(SOURCE)
    lx = bounds["x"][1] - bounds["x"][0]
    ly = bounds["y"][1] - bounds["y"][0]
    xy = tilt[0]
    k2 = phasefield_spectrum(GRID, lx, ly, xy, HEIGHT)

    chosen = None
    attempts = []
    for seed in SEEDS:
        candidate, snapshots = evolve(seed, GRID, k2, max(CHECKPOINTS))
        attempts.append({
            "seed": seed,
            "checkpoints": [{"step": x[3], "phase_field_std": float(np.std(x[4])),
                             "tsc": x[6], "ag": x[7], "accepted": x[0]} for x in snapshots],
        })
        if chosen is None or (candidate[0], candidate[1]) > (chosen[0], chosen[1]):
            chosen = candidate
        if candidate[0]:
            chosen = candidate
            break

    if chosen is None:
        raise RuntimeError("No phase-field candidate was generated")
    accepted, score, seed, step, field, mask, tsc_topology, ag_topology = chosen
    n_tsc_mask = int(mask.sum())
    # Extract a cubic/hexagonal crystal lattice from the checksum-verified planar source.
    c_tsc, d111, ag_z0, tsc_motif, ag_motifs = crystal_lattices(source_atoms, lx, ly, xy)
    tsc_lattice, ag_lattice = [], []
    for repeat in range(math.ceil(HEIGHT / c_tsc)):
        for row in tsc_motif:
            z = float(row[4]) + repeat * c_tsc
            if z < HEIGHT:
                tsc_lattice.append({"type": int(row[1]), "x": float(row[2]), "y": float(row[3]), "z": z})
    for layer in range(math.floor(-ag_z0 / d111) - 1, math.ceil((HEIGHT - ag_z0) / d111) + 1):
        z = ag_z0 + layer * d111
        if 0.0 <= z < HEIGHT:
            for row in ag_motifs[layer % 3]:
                ag_lattice.append({"type": 4, "x": float(row[2]), "y": float(row[3]), "z": z})

    nx, ny, nz = GRID
    mapped = []
    tsc_xyz, ag_xyz = [], []
    for atom in tsc_lattice:
        u, v = frac_xy(atom["x"], atom["y"], lx, ly, xy)
        i, j = min(nx - 1, int(u * nx)), min(ny - 1, int(v * ny))
        k = min(nz - 1, max(0, int(atom["z"] / HEIGHT * nz)))
        if mask[i, j, k]:
            x, y = u * lx + v * xy, v * ly
            mapped.append({**atom, "x": x, "y": y})
            tsc_xyz.append((x, y, atom["z"]))
    for atom in ag_lattice:
        u, v = frac_xy(atom["x"], atom["y"], lx, ly, xy)
        i, j = min(nx - 1, int(u * nx)), min(ny - 1, int(v * ny))
        k = min(nz - 1, max(0, int(atom["z"] / HEIGHT * nz)))
        if not mask[i, j, k]:
            x, y = u * lx + v * xy, v * ly
            mapped.append({**atom, "x": x, "y": y})
            ag_xyz.append((x, y, atom["z"]))

    counts = Counter(a["type"] for a in mapped)
    if not 4000 <= len(mapped) <= 5000:
        raise RuntimeError(f"Mapped structure has {len(mapped)} atoms (target 4000-5000)")
    tsc_xyz = np.asarray(tsc_xyz, float)
    ag_xyz = np.asarray(ag_xyz, float)
    d, _ = min_cross_distances(tsc_xyz, ag_xyz, lx, ly, xy)
    mass_by_type = {1: 47.867, 2: 28.0855, 3: 12.0107, 4: 107.8682}
    total_mass = sum(mass_by_type[a["type"]] for a in mapped)
    ag_mass = counts[4] * mass_by_type[4]
    area = lx * ly
    volume_fraction_tsc_atoms = sum(counts[t] for t in (1, 2, 3)) / len(mapped)
    np.save(ROOT / "phasefield_order_parameter.npy", field.astype(np.float32))
    np.save(ROOT / "phasefield_mask_tsc.npy", mask.astype(np.uint8))
    write_voxel_proxy(ROOT / "AgTi3SiC2_phasefield_voxel_proxy.data", mask, lx, ly, xy)
    map_screen = {
        "status": "direct crystal-site assignment rejected as atomistic input; diagnostic only",
        "method": "TSC and Ag crystal sites classified against the binary phase-field mask",
        "TSC_atoms": int(sum(counts[t] for t in (1, 2, 3))), "Ag_atoms": int(counts[4]),
        "total_atoms": len(mapped), "Ag_mass_fraction_if_counted_as_atoms": ag_mass / total_mass,
        "Ag_TSC_cross_min_distance_A": float(d.min()),
        "Ag_TSC_pairs_below_1p5_A": int((d < 1.5).sum()), "Ag_TSC_pairs_below_2p0_A": int((d < 2.0).sum()),
        "reason_rejected": "minimum cross-phase distance is unphysically short; no pair deletion or relaxation was applied",
        "note": "This diagnostic is not included as a candidate structure and must not be used for molecular dynamics."
    }
    (ROOT / "direct_lattice_mapping_screen.json").write_text(json.dumps(map_screen, indent=2) + "\n")

    metadata = {
        "status": "geometry-only phase-field candidate; no force evaluation, minimization, or dynamics",
        "method": "conserved Cahn-Hilliard phase separation, then quantile assignment of crystalline lattice sites",
        "parameters": {
            "grid": list(GRID), "selected_seed": int(seed), "selected_step": int(step), "dt_arbitrary": DT,
            "mobility_arbitrary": MOBILITY, "kappa_A2_morphology_control": KAPPA_A2,
            "initial_noise_sigma": 0.25, "initial_noise": "random Fourier perturbation restricted to unstable spinodal modes",
            "initial_order_parameter_mean": MEAN_ORDER_PARAMETER,
            "target_TSC_mask_fraction": TARGET_TSC, "mask_voxels_TSC": n_tsc_mask,
            "note": "Mobility, free-energy scale, and gradient coefficient are not calibrated to Ag/Ti3SiC2; only the morphology is used."
        },
        "selected_order_parameter": {"mean": float(field.mean()), "standard_deviation": float(field.std()),
            "minimum": float(field.min()), "maximum": float(field.max())},
        "phase_topology_6_neighbor_periodic_xy_open_z": {"TSC": tsc_topology, "Ag": ag_topology,
            "both_pass_continuity_screen": bool(accepted)},
        "dimensions_A": {"a": [lx, 0.0], "b": [xy, ly], "solid_height": HEIGHT,
            "top_vacuum": VACUUM, "box_z": HEIGHT + VACUUM, "in_plane_angle_deg": 120.0,
            "boundary": "periodic x/y; free top and bottom in z"},
        "crystal_lattices": {"source_sha256": source_hash, "source_file": str(SOURCE),
            "source_is_planar_reference_not_Stage69": True, "TSC_orientation": "(0001) normal +z",
            "Ag_orientation": "(111) normal +z", "TSC_c_repeat_A": c_tsc, "Ag_d111_A": d111,
            "supercells": {"TSC": "12x12", "Ag(111)": "13x13"}},
        "atomistic_mapping_screen": map_screen,
        "voxel_visualization_proxy": {"points": int(mask.size), "TSC_points": int(mask.sum()),
            "Ag_points": int(mask.size - mask.sum()), "not_real_atoms": True,
            "purpose": "load in OVITO to inspect the phase-field morphology only"},
        "attempts": attempts,
        "files": {"voxel_proxy_data": "AgTi3SiC2_phasefield_voxel_proxy.data", "preview": "phasefield_model_preview.png",
            "mask": "phasefield_mask_tsc.npy", "order_parameter": "phasefield_order_parameter.npy",
            "direct_lattice_mapping_screen": "direct_lattice_mapping_screen.json"}
    }
    (ROOT / "audit.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"selected_seed": seed, "selected_step": step, "continuity": accepted,
        "TSC_topology": tsc_topology, "Ag_topology": ag_topology, "atom_count": len(mapped),
        "TSC_atoms": int(sum(counts[t] for t in (1, 2, 3))), "Ag_atoms": int(counts[4]),
        "Ag_mass_fraction": ag_mass / total_mass, "voxel_proxy_points": int(mask.size), "minimum_cross_distance_A": float(d.min()),
        "pairs_below_1p5_A": int((d < 1.5).sum()), "pairs_below_2A": int((d < 2.0).sum())}, indent=2))


if __name__ == "__main__":
    main()
