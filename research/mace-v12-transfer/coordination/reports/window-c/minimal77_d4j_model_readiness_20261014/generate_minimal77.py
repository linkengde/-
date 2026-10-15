#!/usr/bin/env python3
"""Generate the three unrelaxed M77 interface geometries from archived crystallography.

This script is geometry-only. It does not load a calculator or invoke LAMMPS,
DFT, MD, TTM, or a fitting routine.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.geometry import find_mic
from ase.io import read, write


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[5]
DESIGN_DIR = REPO / "research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_ag_interface_minimal_model_design_20261014"
PARENT = REPO / "research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_ag_interface_fit_dft_inputs_20261014/inputs/interface_scan/AGTSC_SI_GAP_2p4A.extxyz"
COD = REPO / "research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif"
COORD_DIR = ROOT / "MINIMAL_MODEL_COORDINATES"

GAPS = (2.16, 2.40, 2.64)
COMMON_A = float(np.sqrt(3.071000001 * 2.889238257))
TSC_A0 = 3.071000001
AG_A0 = 2.889238257
C_Z = 37.9451596
AG_D111 = 2.359
EXPECTED_TSC_Z = np.array([0.0, 2.030772, 3.1534608, 4.434, 5.7145392, 6.837228, 8.868])
EXPECTED_LAYERS = ("Si", "Ti", "C", "Ti", "C", "Ti", "Si")
M = np.array([[2, -1], [1, 3]], dtype=int)  # determinant 7, hexagonal index-7 cell
AREA = 7


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def coset_representatives() -> list[tuple[int, int]]:
    inv = np.linalg.inv(M)
    seen: dict[tuple[int, int], tuple[int, int]] = {}
    for i in range(AREA):
        for j in range(AREA):
            n = np.array([i, j], dtype=float)
            f = np.mod(n @ inv, 1.0)
            key = tuple(np.mod(np.round(f, 10), 1.0))
            seen.setdefault(key, (i, j))
    reps = sorted(seen.values(), key=lambda n: tuple(np.mod(np.round(np.mod(np.array(n) @ inv, 1), 10), 1.0)))
    if len(reps) != AREA:
        raise RuntimeError(f"Expected 7 index cosets, got {len(reps)}")
    return reps


def canonical_layer_offsets(cif: Atoms) -> list[np.ndarray]:
    frac = cif.get_scaled_positions(wrap=True)
    ztargets = (0.25, 0.3645, 0.4278, 0.5, 0.5722, 0.6355, 0.75)
    offsets = []
    for sym, zt in zip(EXPECTED_LAYERS, ztargets):
        candidates = np.where(np.array(cif.get_chemical_symbols()) == sym)[0]
        dz = np.abs(frac[candidates, 2] - zt)
        dz = np.minimum(dz, 1.0 - dz)
        idx = int(candidates[int(np.argmin(dz))])
        if float(np.min(dz)) > 0.01:
            raise RuntimeError(f"Could not identify {sym} layer near z={zt}")
        offsets.append(np.mod(frac[idx, :2], 1.0))
    return offsets


def mic_distance(atoms: Atoms, i: int, j: int) -> float:
    _, d = find_mic(atoms.positions[j] - atoms.positions[i], atoms.cell, atoms.pbc)
    return float(d)


def geometry_rows(atoms: Atoms) -> dict:
    n = len(atoms)
    symbols = np.array(atoms.get_chemical_symbols())
    dm = atoms.get_all_distances(mic=True)
    iu, ju = np.triu_indices(n, 1)
    vals = dm[iu, ju]
    order = np.argsort(vals)
    shortest = []
    for k in order[:20]:
        i, j = int(iu[k]), int(ju[k])
        shortest.append({
            "atom_i": i + 1, "element_i": str(symbols[i]), "region_i": str(atoms.arrays["region"][i]),
            "layer_i": int(atoms.arrays["layer_id"][i]), "atom_j": j + 1,
            "element_j": str(symbols[j]), "region_j": str(atoms.arrays["region"][j]),
            "layer_j": int(atoms.arrays["layer_id"][j]), "distance_A": float(vals[k]),
        })
    pair_min = {}
    for a, b in (("Ag", "Ti"), ("Ag", "Si"), ("Ag", "C")):
        pairs = [(int(i), int(j)) for i in range(n) for j in range(i + 1, n)
                 if {str(symbols[i]), str(symbols[j])} == {a, b}]
        pair_min[f"{a}-{b}"] = min((mic_distance(atoms, i, j), i, j) for i, j in pairs)
    return {
        "min_all_A": float(vals[order[0]]),
        "min_cross_A": {k: {"distance_A": v[0], "atom_i": v[1] + 1, "atom_j": v[2] + 1}
                        for k, v in pair_min.items()},
        "shortest_20": shortest,
        "n_pairs_below_1p5_A": int(np.count_nonzero(vals < 1.5)),
    }


def main() -> None:
    COORD_DIR.mkdir(parents=True, exist_ok=False)
    parent = read(PARENT)
    cif = read(COD)
    parent_tsc = [i for i, region in enumerate(parent.arrays["region"]) if str(region) == "TSC"]
    zsets = []
    for layer in range(1, 8):
        ids = [i for i in parent_tsc if int(parent.arrays["layer_id"][i]) == layer]
        if len(ids) != 19:
            raise RuntimeError(f"Parent layer {layer} has {len(ids)} TSC atoms, expected 19")
        zsets.append(float(np.mean(parent.positions[ids, 2])))
        if Counter(parent[i].symbol for i in ids) != Counter({EXPECTED_LAYERS[layer - 1]: 19}):
            raise RuntimeError(f"Unexpected element in archived parent layer {layer}")
    if not np.allclose(zsets, EXPECTED_TSC_Z, atol=2e-6, rtol=0):
        raise RuntimeError(f"Parent TSC z planes changed: {zsets}")

    # Triangular primitive basis, 120-degree row-vector convention.
    primitive = np.array([[COMMON_A, 0.0], [-COMMON_A / 2.0, np.sqrt(3.0) * COMMON_A / 2.0]])
    super2 = M @ primitive
    cell = np.array([[super2[0, 0], super2[0, 1], 0.0],
                     [super2[1, 0], super2[1, 1], 0.0],
                     [0.0, 0.0, C_Z]])
    invM = np.linalg.inv(M)
    reps = coset_representatives()
    offsets = canonical_layer_offsets(cif)

    base_symbols: list[str] = []
    base_positions: list[list[float]] = []
    regions: list[str] = []
    layers: list[int] = []
    site_ids: list[int] = []
    source_notes: list[str] = []
    base_ids: list[int] = []

    # Substrate order is bottom-to-top planes, then deterministic in-plane sites.
    for layer_idx, (sym, off, z) in enumerate(zip(EXPECTED_LAYERS, offsets, zsets), start=1):
        for site, rep in enumerate(reps, start=1):
            # Keep the integer primitive-cell translation until after mapping to
            # the index-7 supercell; wrapping in the primitive cell would collapse
            # the seven distinct cosets onto one site.
            f_old = np.array(rep, dtype=float) + off
            f_new = np.mod(f_old @ invM, 1.0)
            xy = f_new @ super2
            base_symbols.append(sym)
            base_positions.append([float(xy[0]), float(xy[1]), float(z)])
            regions.append("TSC")
            layers.append(layer_idx)
            site_ids.append(site)
            base_ids.append(len(base_ids) + 1)
            source_notes.append(f"COD9009647_layer{layer_idx}_site{site}")

    # Ag(111) A-B-C-A stacking; A layer is top-Si registered. The layer spacing
    # is kept at the archived parent geometry value, not relaxed for coherency strain.
    stacking = (np.array([0.0, 0.0]), np.array([1 / 3, 2 / 3]),
                np.array([2 / 3, 1 / 3]), np.array([0.0, 0.0]))
    top_si_z = zsets[-1]
    for layer_idx, shift in enumerate(stacking, start=1):
        z = top_si_z + GAPS[1] + (layer_idx - 1) * AG_D111
        for site, rep in enumerate(reps, start=1):
            f_old = np.array(rep, dtype=float) + shift
            f_new = np.mod(f_old @ invM, 1.0)
            xy = f_new @ super2
            base_symbols.append("Ag")
            base_positions.append([float(xy[0]), float(xy[1]), float(z)])
            regions.append("Ag")
            layers.append(layer_idx)
            site_ids.append(site)
            base_ids.append(len(base_ids) + 1)
            source_notes.append(f"Ag111_{'ABCA'[layer_idx-1]}_layer{layer_idx}_site{site}")

    if len(base_symbols) != 77 or Counter(base_symbols) != Counter({"Ag": 28, "Ti": 21, "Si": 14, "C": 14}):
        raise RuntimeError("77-atom composition/order check failed")
    geometry_records = []
    shortest_rows = []
    for gap in GAPS:
        pos = np.asarray(base_positions, float).copy()
        pos[49:, 2] += gap - GAPS[1]
        atoms = Atoms(base_symbols, positions=pos, cell=cell, pbc=(True, True, True))
        atoms.new_array("atom_id", np.arange(1, 78, dtype=int))
        atoms.new_array("region", np.asarray(regions, dtype="U4"))
        atoms.new_array("layer_id", np.asarray(layers, dtype=int))
        atoms.new_array("site_id", np.asarray(site_ids, dtype=int))
        atoms.info.update({
            "structure_id": f"M77_SHARED_STRAIN_GAP_{gap:.2f}A",
            "parent_family": "M77_AREA7_SHARED_STRAIN_DERIVED_FROM_COD9009647_AND_ARCHIVED_217_PARENT",
            "parent_217_sha256": sha(PARENT), "cod_cif_sha256": sha(COD),
            "design_sha256": sha(DESIGN_DIR / "MINIMAL_INTERFACE_DESIGN.md"),
            "gap_definition": "vertical plane-to-plane distance between topmost TSC Si and bottom Ag layer",
            "gap_A": gap, "twist_deg": 0.0, "pbc": "T T T",
            "generator_script_sha256": sha(Path(__file__)),
            "geometry_only_no_relaxation": True,
        })
        tag = f"{gap:.2f}".replace(".", "p")
        path = COORD_DIR / f"Ag_Ti3SiC2_M77_gap_{tag}A.extxyz"
        write(path, atoms, format="extxyz")
        # Fixed-order atom map is shared by all three gap candidates.
        if gap == GAPS[0]:
            with (ROOT / "ATOM_ORDER.csv").open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["atom_id", "element", "region", "layer_id", "site_id", "source_site"])
                for i, row in enumerate(zip(base_symbols, regions, layers, site_ids, source_notes), start=1):
                    w.writerow([i, *row])
        geo = geometry_rows(atoms)
        if geo["n_pairs_below_1p5_A"] != 0:
            raise RuntimeError(f"Unintended sub-1.5 A contact in {atoms.info['structure_id']}")
        if len(atoms) != 77 or Counter(atoms.get_chemical_symbols()) != Counter({"Ag": 28, "Ti": 21, "Si": 14, "C": 14}):
            raise RuntimeError("Post-write geometry identity check failed")
        reread = read(path)
        if len(reread) != 77 or not np.array_equal(reread.arrays["atom_id"], np.arange(1, 78)):
            raise RuntimeError(f"Extxyz round-trip identity failed: {path.name}")
        if not np.array_equal(reread.pbc, np.array([True, True, True])):
            raise RuntimeError(f"Extxyz PBC round-trip failed: {path.name}")
        for rank, item in enumerate(geo["shortest_20"], start=1):
            shortest_rows.append({"structure_id": atoms.info["structure_id"], "rank": rank, **item})
        geometry_records.append({
            "structure_id": atoms.info["structure_id"], "path": str(path.relative_to(ROOT)),
            "sha256": sha(path), "natoms": len(atoms),
            "composition": dict(sorted(Counter(atoms.get_chemical_symbols()).items())),
            "cell_A": cell.tolist(), "pbc": [True, True, True],
            "gap_A": gap, "vacuum_above_top_A": float(C_Z - np.max(atoms.positions[49:, 2])),
            "tSC_plane_z_A": zsets, "Ag_plane_z_A": [float(top_si_z + gap + k * AG_D111) for k in range(4)],
            **geo,
        })

    with (ROOT / "STATIC_GEOMETRY_METRICS.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["structure_id", "path", "sha256", "natoms", "composition", "cell_edge_A",
                    "c_A", "pbc", "gap_A", "vacuum_above_top_A", "min_all_A",
                    "min_Ag_Si_A", "min_Ag_Ti_A", "min_Ag_C_A", "pairs_below_1p5_A",
                    "run0_max_force_eV_A", "run0_status"])
        for s in geometry_records:
            e = s["min_cross_A"]
            w.writerow([s["structure_id"], s["path"], s["sha256"], s["natoms"],
                        json.dumps(s["composition"], sort_keys=True), round(float(np.linalg.norm(cell[0, :2])), 8),
                        C_Z, "T T T", s["gap_A"], s["vacuum_above_top_A"], s["min_all_A"],
                        e["Ag-Si"]["distance_A"], e["Ag-Ti"]["distance_A"], e["Ag-C"]["distance_A"],
                        s["n_pairs_below_1p5_A"], "NOT_MEASURED", "BLOCKED_LAMMPS_NOT_AVAILABLE_AND_CROSS_POTENTIAL_UNVALIDATED"])
    with (ROOT / "SHORTEST_20_DISTANCES.csv").open("w", newline="") as f:
        fields = ["structure_id", "rank", "atom_i", "element_i", "region_i", "layer_i",
                  "atom_j", "element_j", "region_j", "layer_j", "distance_A"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(shortest_rows)

    manifest = {
        "status": "GEOMETRY_GENERATED_UNRELAXED",
        "scope": "coordinates and geometry checks only; no LAMMPS, run 0, MD, TTM, DFT, or fitting",
        "generator_script_sha256": sha(Path(__file__)),
        "design_report_sha256": sha(DESIGN_DIR / "MINIMAL_INTERFACE_DESIGN.md"),
        "COD_CIF_path": str(COD.relative_to(REPO)), "COD_CIF_sha256": sha(COD),
        "archived_217_parent_path": str(PARENT.relative_to(REPO)), "archived_217_parent_sha256": sha(PARENT),
        "natoms": 77, "composition": {"Ag": 28, "Ti": 21, "Si": 14, "C": 14},
        "substrate_plane_order_bottom_to_top": list(EXPECTED_LAYERS),
        "substrate_sites_per_plane": 7,
        "Ag111_layer_stacking_bottom_to_top": ["A", "B", "C", "A"],
        "Ag_sites_per_layer": 7, "registry": "bottom Ag A-layer sites aligned with top Si sites",
        "parent_based_TSC_plane_z_A": zsets, "Ag_interlayer_spacing_A": AG_D111,
        "shared_strain_primitive_a_A": COMMON_A, "index7_matrix": M.tolist(),
        "index7_cell_A": cell.tolist(), "twist_deg": 0.0,
        "TSC_biaxial_strain_pct_from_3p071A": (COMMON_A / TSC_A0 - 1) * 100,
        "Ag_biaxial_strain_pct_from_2p889238A": (COMMON_A / AG_A0 - 1) * 100,
        "pbc": [True, True, True], "gaps_plane_to_plane_A": list(GAPS),
        "c_cell_A": C_Z,
        "limitations": [
            "This coherent M77 compromise is not the archived 12.51983-degree low-mismatch interface.",
            "2.40 A is a geometric reference gap, not an optimized equilibrium separation.",
            "All candidates are unrelaxed and are not DFT labels or force-field validation.",
        ],
        "structures": geometry_records,
    }
    (ROOT / "MODEL_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"structures": geometry_records, "generator_sha256": sha(Path(__file__))}, indent=2))


if __name__ == "__main__":
    main()
