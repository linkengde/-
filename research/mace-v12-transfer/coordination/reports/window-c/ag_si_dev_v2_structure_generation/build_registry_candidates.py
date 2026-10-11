#!/usr/bin/env python3
"""Build four unrelaxed Ag/Ti3SiC2(0001)-Si registry geometries.

This is a geometry-construction and distance-screening script only. It never
creates a calculator or runs DFT, MD, or MACE. The candidate inputs are not
assigned to TRAIN/DEV/TEST and are not claimed to be independent of the
COD-9009647/Kisi parent family.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write


ROOT = Path(__file__).resolve()
while ROOT != ROOT.parent and not (ROOT / "research/mace-v12-transfer").is_dir():
    ROOT = ROOT.parent
if not (ROOT / "research/mace-v12-transfer").is_dir():
    raise RuntimeError("Could not locate repository root")

OUT = Path(__file__).resolve().parent
SOURCE_REL = Path(
    "research/mace-v12-transfer/periodic_interface_v4/"
    "v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif"
)
SOURCE = ROOT / SOURCE_REL
EXPECTED_SOURCE_SHA256 = "6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3"

# Source-backed, unrelaxed screening geometry choices. See REGISTRY_DESIGN.md.
MAX_COD_A = 17.6235
MAX_A_INTERFACE_A = 3.071
AG_FCC_A_A = 4.09
AG_111_NN_A = AG_FCC_A_A / np.sqrt(2.0)
# Keep the constructed Ag(111) slab geometrically fcc under the coherent in-
# plane match: scale the ideal cubic Ag lattice isotropically to sqrt(2)*aMAX.
AG_MATCHED_FCC_A_A = np.sqrt(2.0) * MAX_A_INTERFACE_A
AG_111_D_A = AG_MATCHED_FCC_A_A / np.sqrt(3.0)
SUBSTRATE_BLOCKS = 2
SUBSTRATE_THICKNESS_A = SUBSTRATE_BLOCKS * MAX_COD_A / 2.0
PLANE_GAP_A = 2.60
VACUUM_BOTTOM_A = 18.0
VACUUM_TOP_A = 18.0
AG_LAYERS = 7
AG_STACK = [(0.0, 0.0), (2 / 3, 1 / 3), (1 / 3, 2 / 3)]

REGISTRIES = [
    {
        "id": "R1",
        "filename": "Ag_Ti3SiC2_Si_R1.xyz",
        "name": "top / OT-like",
        "frac": (0.0, 0.0),
        "definition": "bottom Ag site directly above the terminating Si site",
    },
    {
        "id": "R2",
        "filename": "Ag_Ti3SiC2_Si_R2.xyz",
        "name": "bridge / MT-like",
        "frac": (0.5, 0.0),
        "definition": "bottom Ag site at the midpoint of adjacent terminating-Si lattice sites along a1",
    },
    {
        "id": "R3",
        "filename": "Ag_Ti3SiC2_Si_R3.xyz",
        "name": "threefold hollow / HCP-like",
        "frac": (2 / 3, 1 / 3),
        "definition": "hollow translation placing bottom Ag above the second-layer Ti lateral site",
    },
    {
        "id": "R4",
        "filename": "Ag_Ti3SiC2_Si_R4.xyz",
        "name": "alternate threefold hollow / FCC-like",
        "frac": (1 / 3, 2 / 3),
        "definition": "opposite hollow translation, above the third-layer C lateral site",
    },
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def build_substrate():
    if sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise RuntimeError("COD source CIF SHA-256 mismatch")
    bulk = read(SOURCE, format="cif")
    if len(bulk) != 12 or bulk.get_chemical_formula() != "C4Si2Ti6":
        raise RuntimeError("Unexpected COD parent composition or unit-cell multiplicity")
    cell0 = np.asarray(bulk.cell)
    if not np.allclose(cell0[2], [0.0, 0.0, MAX_COD_A], atol=1e-5):
        raise RuntimeError("Unexpected parent c parameter")

    frac = bulk.get_scaled_positions(wrap=True)
    syms = bulk.get_chemical_symbols()
    rows = []
    z_lo_frac = 0.25
    # COD's 12-atom conventional cell has two Si-to-Si Ti3C2 blocks per c.
    z_hi_frac = z_lo_frac + SUBSTRATE_BLOCKS / 2.0
    for rep_c in range(SUBSTRATE_BLOCKS + 1):
        for source_i, (sym, xyzf) in enumerate(zip(syms, frac), start=1):
            zf = float(xyzf[2]) + rep_c
            if z_lo_frac - 1e-7 <= zf <= z_hi_frac + 1e-7:
                z = (zf - z_lo_frac) * MAX_COD_A
                rows.append({
                    "symbol": sym,
                    "fx": float(xyzf[0] % 1.0),
                    "fy": float(xyzf[1] % 1.0),
                    "z": z,
                    "source_i": source_i,
                    "source_repeat_c": rep_c,
                })
    rows.sort(key=lambda r: (round(r["z"], 7), round(r["fx"], 7), round(r["fy"], 7), r["symbol"]))
    if len(rows) != 13:
        raise RuntimeError(f"Expected 13 substrate planes/atoms per 1x1 cell, found {len(rows)}")
    if [r["symbol"] for r in rows] != ["Si", "Ti", "C", "Ti", "C", "Ti", "Si", "Ti", "C", "Ti", "C", "Ti", "Si"]:
        raise RuntimeError("Unexpected Si-terminated layer sequence")
    if not np.isclose(rows[0]["z"], 0.0, atol=1e-7) or not np.isclose(rows[-1]["z"], SUBSTRATE_THICKNESS_A, atol=1e-5):
        raise RuntimeError("Si termination planes do not bound the selected substrate slab")
    return rows, sha256(SOURCE)


def build_atoms(registry, substrate_rows, source_hash):
    a = MAX_A_INTERFACE_A
    a1 = np.array([a, 0.0, 0.0])
    a2 = np.array([-a / 2.0, np.sqrt(3.0) * a / 2.0, 0.0])
    cell_z = VACUUM_BOTTOM_A + SUBSTRATE_THICKNESS_A + PLANE_GAP_A + (AG_LAYERS - 1) * AG_111_D_A + VACUUM_TOP_A
    cell = np.array([a1, a2, [0.0, 0.0, cell_z]])

    syms, positions = [], []
    atom_ids, subsystem, layer_index, source_atom_index, source_repeat_c = [], [], [], [], []
    for plane_i, row in enumerate(substrate_rows):
        syms.append(row["symbol"])
        positions.append(row["fx"] * a1 + row["fy"] * a2 + np.array([0.0, 0.0, VACUUM_BOTTOM_A + row["z"]]))
        atom_ids.append(len(atom_ids) + 1)
        subsystem.append(0)
        layer_index.append(plane_i + 1)
        source_atom_index.append(row["source_i"])
        source_repeat_c.append(row["source_repeat_c"])

    t = np.array(registry["frac"], dtype=float)
    z_interface = VACUUM_BOTTOM_A + SUBSTRATE_THICKNESS_A + PLANE_GAP_A
    for layer in range(AG_LAYERS):
        stack = np.array(AG_STACK[layer % 3])
        xyf = (t + stack) % 1.0
        syms.append("Ag")
        positions.append(xyf[0] * a1 + xyf[1] * a2 + np.array([0.0, 0.0, z_interface + layer * AG_111_D_A]))
        atom_ids.append(len(atom_ids) + 1)
        subsystem.append(1)
        layer_index.append(layer + 1)
        source_atom_index.append(0)
        source_repeat_c.append(-1)

    atoms = Atoms(symbols=syms, positions=positions, cell=cell, pbc=[True, True, False])
    atoms.new_array("atom_id", np.array(atom_ids, dtype=int))
    atoms.new_array("subsystem", np.array(subsystem, dtype=int))  # 0=MAX substrate, 1=Ag
    atoms.new_array("layer_index", np.array(layer_index, dtype=int))
    atoms.new_array("source_atom_index", np.array(source_atom_index, dtype=int))
    atoms.new_array("source_repeat_c", np.array(source_repeat_c, dtype=int))
    strain_pct = 100.0 * (MAX_A_INTERFACE_A / AG_111_NN_A - 1.0)
    atoms.info.update({
        "structure_id": registry["filename"].removesuffix(".xyz"),
        "config_type": "unrelaxed_registry_geometry_candidate",
        "source_family": "COD-9009647_Kisi-1998-derived; same crystal parent family as prior scored Ag-Si structures",
        "structure_role": "GEOMETRY_ONLY_NOT_ASSIGNED_TRAIN_DEV_TEST",
        "registry": registry["name"],
        "registry_fractional_xy": f"{registry['frac'][0]:.8f},{registry['frac'][1]:.8f}",
        "mean_plane_gap_A": float(PLANE_GAP_A),
        "ag_lattice_constant_fcc_A": float(AG_FCC_A_A),
        "ag111_unstrained_nn_A": float(AG_111_NN_A),
        "ag111_coherently_strained_fcc_a_A": float(AG_MATCHED_FCC_A_A),
        "ag111_coherently_strained_interlayer_spacing_A": float(AG_111_D_A),
        "ag_inplane_coherent_strain_percent": float(strain_pct),
        "substrate_cod_parent_sha256": source_hash,
        "unrelaxed_geometry": True,
        "dft_authorized": False,
    })
    return atoms


def min_pair_distance(atoms, source_indices, target_indices):
    return min(atoms.get_distance(int(i), int(j), mic=True) for i in source_indices for j in target_indices)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    substrate_rows, source_hash = build_substrate()
    substrate_formula = " ".join(f"{s}:{sum(r['symbol']==s for r in substrate_rows)}" for s in ("Ti", "Si", "C"))
    strain_pct = 100.0 * (MAX_A_INTERFACE_A / AG_111_NN_A - 1.0)
    results = []
    for reg in REGISTRIES:
        atoms = build_atoms(reg, substrate_rows, source_hash)
        path = OUT / reg["filename"]
        write(path, atoms, format="extxyz")
        check = read(path, format="extxyz")
        if len(check) != 20 or check.get_chemical_formula() != "C4Ag7Si3Ti6":
            raise RuntimeError(f"Unexpected atom count/formula in {path.name}")
        if list(check.pbc) != [True, True, False]:
            raise RuntimeError(f"Unexpected PBC in {path.name}")
        if not np.allclose(check.cell, atoms.cell, atol=1e-10):
            raise RuntimeError(f"Cell changed during XYZ roundtrip for {path.name}")

        ag = [i for i, s in enumerate(check.get_chemical_symbols()) if s == "Ag"]
        si = [i for i, s in enumerate(check.get_chemical_symbols()) if s == "Si"]
        ti = [i for i, s in enumerate(check.get_chemical_symbols()) if s == "Ti"]
        carbon = [i for i, s in enumerate(check.get_chemical_symbols()) if s == "C"]
        d_ag_si = min_pair_distance(check, ag, si)
        d_ag_ti = min_pair_distance(check, ag, ti)
        d_ag_c = min_pair_distance(check, ag, carbon)
        severe = min(d_ag_si, d_ag_ti, d_ag_c) < 2.0 - 1e-10
        results.append({
            "registry_id": reg["id"],
            "structure_file": reg["filename"],
            "registry_name": reg["name"],
            "definition": reg["definition"],
            "Ag_coverage": "1 Ag atom per 1x1 coherent cell per layer; full 1x1 coverage",
            "Ag_bottom_layer_fractional_shift": f"{reg['frac'][0]:.8f},{reg['frac'][1]:.8f}",
            "n_atoms": len(check),
            "formula": check.get_chemical_formula(),
            "a_substrate_cell_A": f"{MAX_A_INTERFACE_A:.6f}",
            "Ag_fcc_a_source_A": f"{AG_FCC_A_A:.6f}",
            "Ag111_matched_fcc_a_A": f"{AG_MATCHED_FCC_A_A:.6f}",
            "Ag111_interlayer_spacing_A": f"{AG_111_D_A:.6f}",
            "mean_plane_gap_A": f"{PLANE_GAP_A:.6f}",
            "min_Ag_Si_A": f"{d_ag_si:.6f}",
            "min_Ag_Ti_A": f"{d_ag_ti:.6f}",
            "min_Ag_C_A": f"{d_ag_c:.6f}",
            "Ag_Si_cutoff_pass_ge_2A": str(d_ag_si >= 2.0).upper(),
            "Ag_Ti_C_severe_overlap_pass_ge_2A": str(min(d_ag_ti, d_ag_c) >= 2.0).upper(),
            "all_pair_geometry_screen": "FAIL" if severe else "PASS",
            "Ag_inplane_strain_pct_tensile": f"{strain_pct:.6f}",
            "source_family_status": "COD-9009647-derived; not an independent DEV family",
            "dft_status": "NOT_RUN_NOT_AUTHORIZED",
            "sha256": sha256(path),
        })
        if severe:
            raise RuntimeError(f"Severe Ag/substrate overlap in {path.name}")

    with (OUT / "REGISTRY_SCREENING.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(results)

    manifest = {
        "status": "READY_FOR_DFT_SCREENING",
        "dft_execution_authorized": False,
        "data_role": "GEOMETRY_ONLY_NOT_ASSIGNED_TRAIN_DEV_TEST",
        "independent_dev_family": False,
        "independence_note": "Substrate derives from COD 9009647/Kisi 1998, the same crystallographic parent family as prior scored Ag-Si geometries; do not treat these as independent DEV evidence.",
        "generator_sha256": sha256(Path(__file__).resolve()),
        "source_cif": str(SOURCE_REL),
        "source_cif_sha256": source_hash,
        "geometry_parameters": {
            "max_interface_inplane_a_A": MAX_A_INTERFACE_A,
            "substrate_thickness_A": SUBSTRATE_THICKNESS_A,
            "substrate_atomic_planes": len(substrate_rows),
            "ag_fcc_source_a_A": AG_FCC_A_A,
            "ag111_unstrained_nn_A": AG_111_NN_A,
            "ag111_coherently_strained_fcc_a_A": AG_MATCHED_FCC_A_A,
            "ag111_interlayer_spacing_A": AG_111_D_A,
            "ag_layers": AG_LAYERS,
            "mean_plane_gap_A": PLANE_GAP_A,
            "vacuum_bottom_A": VACUUM_BOTTOM_A,
            "vacuum_top_A": VACUUM_TOP_A,
            "pbc": [True, True, False],
        },
        "substrate_composition_per_candidate": substrate_formula,
        "candidates": results,
    }
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
