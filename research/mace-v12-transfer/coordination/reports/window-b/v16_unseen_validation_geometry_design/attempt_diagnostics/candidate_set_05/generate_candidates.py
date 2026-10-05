#!/usr/bin/env python3
"""Geometry-only design and deduplication of three proposed V16 holdout inputs.

Reads dataset geometries, reserved input geometries, prior label-free proposals and
the historical V15 localization report. It never reads calculation output folders,
DFT reference forces/energies, or loads a model.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

PROJECT = Path("research/mace-v12-transfer/periodic_interface_v4")
WINDOW_B = Path("research/mace-v12-transfer/coordination/reports/window-b")
REPORT = WINDOW_B / "v16_unseen_validation_geometry_design"
INTERFACE = {
    "AgC": "AgC_registry_strain_acq_v14_01_periodic_PW_PBE_energy_force_v14_train",
    "AgSi": "AgSi_registry_strain_acq_v14_01_periodic_PW_PBE_energy_force_v14_train",
    "AgTi": "AgTi_registry_probe_v13_01_periodic_PW_PBE_energy_force_v14_train",
}
PAIR_FLOOR_MARGIN_A = 0.05
NEAR_REGISTRY_RMS_A = 0.15
NEAR_FRAMEWORK_RMS_A = 0.15


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def id_map(atoms):
    ids = np.asarray(atoms.arrays.get("lammps_id", []), dtype=int)
    if len(ids) != len(atoms) or len(set(ids.tolist())) != len(atoms):
        return None
    return {int(k): i for i, k in enumerate(ids)}


def formula(a):
    return a.get_chemical_formula()


def pbc_cell_match(a, b):
    return np.array_equal(a.pbc, b.pbc) and np.allclose(a.cell.array, b.cell.array, atol=1e-6, rtol=0)


def geometry_compare(a, b):
    """Translation-aligned geometry distance, ID first and species assignment fallback."""
    from ase.geometry import find_mic

    if len(a) != len(b) or formula(a) != formula(b) or not pbc_cell_match(a, b):
        return None
    am, bm = id_map(a), id_map(b)
    aligned_by_id = am is not None and bm is not None and set(am) == set(bm)
    if aligned_by_id:
        order = sorted(am)
        di = np.asarray([a.positions[am[k]] - b.positions[bm[k]] for k in order])
        symbols = np.asarray([a[am[k]].symbol for k in order])
        di, _ = find_mic(di, a.cell, pbc=a.pbc)
    else:
        # Conservative symmetry-aware fallback for same-composition files without shared IDs.
        from scipy.optimize import linear_sum_assignment
        need(sorted(a.get_chemical_symbols()) == sorted(b.get_chemical_symbols()), "same formula but element multiset differs")
        # Remove a common translation before species-wise minimum-cost assignment.
        aa = a.positions - np.mean(a.positions, axis=0)
        bb = b.positions - np.mean(b.positions, axis=0)
        pairs = []
        for element in sorted(set(a.get_chemical_symbols())):
            ai = np.asarray([i for i, x in enumerate(a.get_chemical_symbols()) if x == element], dtype=int)
            bi = np.asarray([i for i, x in enumerate(b.get_chemical_symbols()) if x == element], dtype=int)
            diff = aa[ai, None, :] - bb[None, bi, :]
            wrapped, d = find_mic(diff.reshape(-1, 3), a.cell, pbc=a.pbc)
            cost = np.asarray(d).reshape(len(ai), len(bi))
            ri, ci = linear_sum_assignment(cost)
            pairs.extend((int(i), int(j)) for i, j in zip(ai[ri], bi[ci]))
        di = np.asarray([a.positions[i] - b.positions[j] for i, j in pairs])
        symbols = np.asarray([a[i].symbol for i, _ in pairs])
        di, _ = find_mic(di, a.cell, pbc=a.pbc)
    fw = symbols != "Ag"
    common = np.median(di[fw], axis=0) if fw.any() else np.median(di, axis=0)
    di = di - common
    mag = np.linalg.norm(di, axis=1)
    ag = symbols == "Ag"
    return {
        "method": "persistent_lammps_id" if aligned_by_id else "species_assignment_fallback",
        "all_atom_RMS_A": float(np.sqrt(np.mean(mag**2))),
        "Ag_registry_RMS_A": float(np.sqrt(np.mean(mag[ag]**2))) if ag.any() else 0.0,
        "framework_RMS_A": float(np.sqrt(np.mean(mag[fw]**2))) if fw.any() else 0.0,
        "max_atom_displacement_A": float(mag.max()),
        "exact_duplicate": bool(mag.max() <= 1e-5),
        "near_duplicate": bool((mag.max() <= 1e-5) or (np.sqrt(np.mean(mag[ag]**2)) <= NEAR_REGISTRY_RMS_A and np.sqrt(np.mean(mag[fw]**2)) <= NEAR_FRAMEWORK_RMS_A)),
    }


def pair_minima(atoms):
    d = atoms.get_all_distances(mic=True)
    symbols = atoms.get_chemical_symbols()
    out = {}
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            pair = "-".join(sorted((symbols[i], symbols[j])))
            value = float(d[i, j])
            if pair not in out or value < out[pair]["distance_A"]:
                ids = np.asarray(atoms.arrays.get("lammps_id", np.arange(len(atoms))), dtype=int)
                out[pair] = {"distance_A": value, "ids": [int(ids[i]), int(ids[j])], "elements": [symbols[i], symbols[j]]}
    return out


def central_pair(atoms):
    need("central_pair" in atoms.arrays and "lammps_id" in atoms.arrays, "source lacks central-pair marker/IDs")
    marked = np.flatnonzero(np.asarray(atoms.arrays["central_pair"], dtype=bool))
    need(len(marked) == 2, "expected two marked atoms")
    ag = [int(i) for i in marked if atoms[i].symbol == "Ag"]
    x = [int(i) for i in marked if atoms[i].symbol != "Ag"]
    need(len(ag) == len(x) == 1, "marked pair must be Ag-X")
    return ag[0], x[0]


def clean_candidate(source, label, translation, source_path, source_sha, geometry_sha):
    a = source.copy()
    ag_indices = [i for i, atom in enumerate(a) if atom.symbol == "Ag"]
    a.positions[ag_indices] += np.asarray(translation, dtype=float)
    allowed_arrays = {"lammps_id", "lammps_type", "central_pair"}
    for key in list(a.arrays):
        if key not in allowed_arrays and key not in {"numbers", "positions"}:
            del a.arrays[key]
    a.info = {
        "config_type": label,
        "source_config_type": str(source.info.get("config_type", "")),
        "source_dataset_path": source_path,
        "source_dataset_sha256": source_sha,
        "source_frame_geometry_sha256": geometry_sha,
        "registry_translation_A": np.asarray(translation).tolist(),
        "candidate_role": "unassigned_withheld_validation_proposal",
        "proposal_only_no_DFT_label": True,
    }
    a.calc = None
    return a


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()
    repo, out = args.repo_root.resolve(), args.output_dir.resolve()
    need(out.is_relative_to(repo / REPORT) and out != repo / REPORT, "output must be a child of this task report")
    need(not out.exists() or (out.is_dir() and not any(out.iterdir())), "nonempty output refused")
    out.mkdir(parents=True, exist_ok=True)
    evidence_hashes, inventory, screening = {}, [], []

    def record(path, role):
        path = Path(path).resolve()
        need(path.is_relative_to(repo) and path.is_file(), f"missing/outside-repository source: {path}")
        relative = str(path.relative_to(repo))
        digest = sha(path)
        evidence_hashes[relative] = digest
        return path, digest

    try:
        from ase.io import read, write
        from ase.geometry import find_mic

        # Dataset split geometries; labels in existing archives are never accessed.
        data_files = sorted((repo / PROJECT).glob("mace_periodic_v*/data/*.extxyz"))
        split_geometries, all_known, official_input_geometries = [], [], []
        for path in data_files:
            if path.name not in {"train.extxyz", "valid.extxyz", "test.extxyz"}:
                continue
            p, digest = record(path, "current_split_geometry")
            frames = read(p, index=":")
            inventory.append({"path": str(p.relative_to(repo)), "sha256": digest, "frames": len(frames), "role": "existing train/valid/test geometry only"})
            for index, frame in enumerate(frames):
                frame.calc = None
                split_geometries.append((str(p.relative_to(repo)), index, frame))
                all_known.append((str(p.relative_to(repo)), index, frame))

        # Frozen/assigned DFT input manifests only. Never traverse any calculation directory.
        input_dirs = [
            PROJECT / "pbe_interface_v13_holdouts",
            PROJECT / "pbe_interface_v14_main_holdouts",
            PROJECT / "pbe_interface_v14_parallel_acquisition",
            PROJECT / "pbe_interface_v15_registry_holdouts",
            PROJECT / "pbe_interface_v15_targeted_acquisition",
            PROJECT / "pbe_interface_v16_parallel_acquisition",
            PROJECT / "pbe_interface_v16_main_acquisition",
        ]
        planned_acquisitions = []
        for rel in input_dirs:
            manifest_path = repo / rel / "input_manifest.json"
            if not manifest_path.is_file():
                continue
            p, digest = record(manifest_path, "input_manifest")
            manifest = json.loads(p.read_text())
            inventory.append({"path": str(p.relative_to(repo)), "sha256": digest, "frames": len(manifest.get("records", [])), "role": "input-manifest-only"})
            for rec in manifest.get("records", []):
                ipath = (repo / rel / rec["input"]).resolve()
                fp, fh = record(ipath, "input_geometry")
                need(fh == rec["input_sha256"], f"frozen input hash mismatch: {rec['label']}")
                frame = read(fp)
                need("REF_energy" not in frame.info and "PW_PBE_energy_eV" not in frame.info and not any(k in frame.arrays for k in ("REF_forces", "PW_PBE_forces")), f"input contains labels: {rec['label']}")
                frame.calc = None
                item = (str(fp.relative_to(repo)), 0, frame)
                all_known.append(item)
                official_input_geometries.append(item)
                if rel.name in {"pbe_interface_v16_parallel_acquisition", "pbe_interface_v16_main_acquisition"}:
                    planned_acquisitions.append({"label": rec["label"], "sha256": fh, "path": str(fp.relative_to(repo)), "owner": rec.get("proposed_owner"), "candidate_id": rec.get("candidate_id"), "atoms": len(frame), "formula": formula(frame)})

        # Earlier B label-free geometry proposals also count as prior designs.
        proposal_dirs = [
            repo / WINDOW_B / "v15_registry_candidates",
            repo / WINDOW_B / "v15_wide_span_sampling_design/candidates",
        ]
        for directory in proposal_dirs:
            for fp in sorted(directory.glob("*.extxyz")):
                p, digest = record(fp, "prior_label_free_candidate")
                frames = read(p, index=":")
                inventory.append({"path": str(p.relative_to(repo)), "sha256": digest, "frames": len(frames), "role": "prior label-free geometry proposal"})
                for index, frame in enumerate(frames):
                    frame.calc = None
                    all_known.append((str(p.relative_to(repo)), index, frame))

        local_path = repo / WINDOW_B / "v15_force_error_localization/analysis_01/localization.json"
        lp, lhash = record(local_path, "historical_force_localization")
        localization = json.loads(lp.read_text())
        historical_targets = {}
        for source in localization["sources"]:
            label = source["label"].split("_")[0]
            # Prioritize non-Ag local atoms, not only the marked-pair scalar projection.
            historical_targets[label] = [{"id": int(row["id"]), "element": row["element"], "error_norm_eV_A": float(row["force_error_norm_eV_A"])} for row in source["top_10_atom_errors"] if row["element"] != "Ag"][:4]

        # Source mapping and geometry hashes from the earlier label-free screen.
        parent_screen_path = repo / WINDOW_B / "v15_registry_candidates/candidate_screening_manifest.json"
        pp, phash = record(parent_screen_path, "parent_lineage_geometry_manifest")
        parent_screen = json.loads(pp.read_text())
        parent_meta = parent_screen["source_dataset"]["parent_configs"]
        parent_dataset_path = repo / PROJECT / "mace_periodic_v14_interface_energy/data/train.extxyz"
        parent_file, parent_sha = record(parent_dataset_path, "parent_dataset_geometry")
        parent_frames = read(parent_file, index=":")
        parents = {}
        for kind, cfg in INTERFACE.items():
            rows = [(i, a) for i, a in enumerate(parent_frames) if a.info.get("config_type") == cfg]
            need(len(rows) == 1, f"parent frame not unique: {kind}")
            parents[kind] = rows[0]

        # Include existing V15 holdout-input and planned V16 registry translation vectors.
        existing_registry_vectors = {kind: [] for kind in INTERFACE}
        for path_text, index, frame in all_known:
            info = frame.info
            cfg = str(info.get("source_config_type", ""))
            for kind, source_cfg in INTERFACE.items():
                if cfg == source_cfg and "registry_translation_A" in info:
                    existing_registry_vectors[kind].append({"path": path_text, "translation_A": np.asarray(info["registry_translation_A"], dtype=float).tolist()})
        for kind in INTERFACE:
            rec = parent_meta[kind]
            for input_rec in planned_acquisitions:
                if input_rec["label"].startswith(kind) and "registry" in input_rec["label"]:
                    frame = read(repo / input_rec["path"])
                    ag = np.asarray([i for i, a in enumerate(frame) if a.symbol == "Ag"], dtype=int)
                    source = parents[kind][1]
                    if set(id_map(source) or {}) == set(id_map(frame) or {}):
                        sm, fm = id_map(source), id_map(frame)
                        deltas, _ = find_mic(np.asarray([frame.positions[fm[k]] - source.positions[sm[k]] for k in sorted(sm) if source[sm[k]].symbol == "Ag"]), source.cell, pbc=source.pbc)
                        shift = np.median(deltas, axis=0)
                        existing_registry_vectors[kind].append({"path": input_rec["path"], "translation_A": shift.tolist()})

        # Species-aware contact floors from the same formula's current training/validation/test
        # and reserved/acquisition inputs. These are empirical guardrails, not bond definitions.
        pools = {kind: [x for x in all_known if formula(x[2]) == formula(parents[kind][1])] for kind in INTERFACE}
        chemical_pools = {kind: [x for x in (split_geometries + official_input_geometries) if formula(x[2]) == formula(parents[kind][1])] for kind in INTERFACE}
        observed_minima, central_ranges = {}, {}
        for kind in INTERFACE:
            by_pair = {}
            marked_values = []
            for path_text, index, frame in chemical_pools[kind]:
                for pair, value in pair_minima(frame).items():
                    if pair not in by_pair or value["distance_A"] < by_pair[pair]["distance_A"]:
                        by_pair[pair] = {**value, "path": path_text, "frame": index}
                try:
                    ia, ix = central_pair(frame)
                    marked_values.append({"distance_A": float(frame.get_distance(ia, ix, mic=True)), "path": path_text, "frame": index})
                except Exception:
                    pass
            observed_minima[kind] = by_pair
            central_ranges[kind] = [{"distance_A": min(marked_values, key=lambda x: x["distance_A"])["distance_A"], "source": min(marked_values, key=lambda x: x["distance_A"])}, {"distance_A": max(marked_values, key=lambda x: x["distance_A"])["distance_A"], "source": max(marked_values, key=lambda x: x["distance_A"])}] if marked_values else None

        # Search a modest 2D grid for a distinct lateral Ag registry. The framework
        # remains fixed; thus its force response samples a new Ag-neighbor environment.
        magnitudes = [0.35, 0.45, 0.55, 0.65, 0.75, 0.85]
        angles = list(range(0, 360, 15))
        candidates = []
        for kind, cfg in INTERFACE.items():
            frame_index, source = parents[kind]
            meta = parent_meta[kind]
            ia, ix = central_pair(source)
            pair_ids = np.asarray(source.arrays["lammps_id"], dtype=int)
            ag_indices = [i for i, atom in enumerate(source) if atom.symbol == "Ag"]
            # Established interface normal is z; search the xy registry plane.
            axis = source.positions[ix] - source.positions[ia]
            lateral_axis = np.cross(axis, np.asarray([0.0, 0.0, 1.0]))
            lateral_axis[2] = 0.0
            need(np.linalg.norm(lateral_axis) > 1e-9, f"cannot define lateral axis: {kind}")
            lateral_axis = lateral_axis / np.linalg.norm(lateral_axis)
            targets = [t for t in historical_targets.get(kind, []) if t["id"] in set(pair_ids.tolist())]
            best = None
            candidate_checks = []
            for magnitude in magnitudes:
                for angle in angles:
                    theta = math.radians(angle)
                    rotated = np.asarray([math.cos(theta) * lateral_axis[0] - math.sin(theta) * lateral_axis[1], math.sin(theta) * lateral_axis[0] + math.cos(theta) * lateral_axis[1], 0.0])
                    translation = magnitude * rotated
                    label = f"Ag{kind[2:]}_v16_unseen_validation_01"
                    candidate = clean_candidate(source, label, translation, str(parent_dataset_path.relative_to(repo)), parent_sha, meta["source_frame_geometry_sha256"])
                    minima = pair_minima(candidate)
                    floor_fail = []
                    for pair, floor_record in observed_minima[kind].items():
                        if pair not in minima:
                            continue
                        floor = floor_record["distance_A"] - PAIR_FLOOR_MARGIN_A
                        if minima[pair]["distance_A"] < floor - 1e-10:
                            floor_fail.append({"pair": pair, "candidate_min_A": minima[pair]["distance_A"], "empirical_min_A": floor_record["distance_A"], "screen_floor_A": floor, "empirical_minimum_source": floor_record})
                    central_d = float(candidate.get_distance(ia, ix, mic=True))
                    crange = central_ranges[kind]
                    gap_fail = crange is None or central_d < crange[0]["distance_A"] - 0.20 or central_d > crange[1]["distance_A"] + 0.20
                    comparisons = []
                    exact = []
                    near = []
                    for path_text, idx, old in pools[kind]:
                        comparison = geometry_compare(candidate, old)
                        if comparison is None:
                            continue
                        item = {"path": path_text, "frame": idx, **comparison}
                        comparisons.append(item)
                        if comparison["exact_duplicate"]:
                            exact.append(item)
                        if comparison["near_duplicate"]:
                            near.append(item)
                    # Candidate must differ from every same-composition known frame under
                    # the original registry-RMS near-duplicate rule, plus exact check.
                    duplicate_fail = bool(exact or near)
                    target_changes = []
                    sm = id_map(source); cm = id_map(candidate)
                    for target in targets:
                        k = target["id"]
                        i = cm[k]
                        old_min = min(float(source.get_distance(i, j, mic=True)) for j in ag_indices if j != i) if source[i].symbol != "Ag" else 0.0
                        new_min = min(float(candidate.get_distance(i, j, mic=True)) for j in ag_indices) if candidate[i].symbol != "Ag" else 0.0
                        target_changes.append({"id": k, "element": target["element"], "historical_error_norm_eV_A": target["error_norm_eV_A"], "nearest_Ag_distance_parent_A": old_min, "nearest_Ag_distance_candidate_A": new_min, "nearest_Ag_distance_change_A": new_min - old_min})
                    target_effect = max((abs(x["nearest_Ag_distance_change_A"]) for x in target_changes), default=0.0)
                    nearest = min(comparisons, key=lambda x: x["all_atom_RMS_A"]) if comparisons else None
                    novelty = float(nearest["all_atom_RMS_A"]) if nearest else 1.0
                    valid = not floor_fail and not gap_fail and not duplicate_fail
                    check = {"magnitude_A": magnitude, "angle_from_established_lateral_axis_deg": angle, "translation_A": translation.tolist(), "marked_pair_distance_A": central_d, "species_pair_minima": minima, "empirical_species_floor_failures": floor_fail, "marked_pair_range_fail": gap_fail, "exact_or_near_duplicate_fail": duplicate_fail, "exact_hits": exact, "near_hits": near, "nearest_known_geometry": nearest, "target_neighbor_changes": target_changes, "max_target_Ag_distance_change_A": target_effect, "valid": valid, "selection_score": novelty + 0.25 * target_effect}
                    candidate_checks.append(check)
                    if valid and (best is None or check["selection_score"] > best[0]):
                        best = (check["selection_score"], candidate, check)
            need(best is not None, f"no candidate passes geometry, species-aware contact and split-dedup screens for {kind}")
            score, selected, chosen = best
            # Recheck exact atom ordering, lineage, PBC, and absent labels.
            need(np.array_equal(selected.numbers, source.numbers) and np.array_equal(selected.arrays["lammps_id"], source.arrays["lammps_id"]), f"candidate atom order/IDs changed: {kind}")
            need(np.array_equal(selected.pbc, source.pbc) and np.allclose(selected.cell.array, source.cell.array, atol=1e-12, rtol=0), f"candidate cell/PBC changed: {kind}")
            need("REF_energy" not in selected.info and not any(k in selected.arrays for k in ("REF_forces", "PW_PBE_forces", "forces")), f"candidate contains labels: {kind}")
            outpath = out / f"{kind}_v16_unseen_validation_01.extxyz"
            write(outpath, selected, format="extxyz")
            # Read back and verify that writing preserved geometry/IDs while not adding labels.
            check = read(outpath)
            need(np.array_equal(check.numbers, selected.numbers) and np.array_equal(check.arrays["lammps_id"], selected.arrays["lammps_id"]), f"serialized candidate atom identity mismatch: {kind}")
            # ASE extxyz writes Cartesian coordinates to eight decimal places.
            need(np.allclose(check.positions, selected.positions, atol=1e-8, rtol=0) and np.allclose(check.cell.array, selected.cell.array, atol=1e-12, rtol=0) and np.array_equal(check.pbc, selected.pbc), f"serialized candidate geometry changed beyond extxyz precision: {kind}")
            need("REF_energy" not in check.info and not any(k in check.arrays for k in ("REF_forces", "PW_PBE_forces", "forces")), f"serialized candidate contains labels: {kind}")
            target_neighbor_before_after = chosen["target_neighbor_changes"]
            species_contact_audit = {}
            for pair, minimum in chosen["species_pair_minima"].items():
                observed = observed_minima[kind].get(pair)
                floor = None if observed is None else float(observed["distance_A"] - PAIR_FLOOR_MARGIN_A)
                species_contact_audit[pair] = {"candidate_minimum_A": minimum["distance_A"], "smallest_observed_same_composition_A": None if observed is None else observed["distance_A"], "empirical_minimum_source": observed, "empirical_screen_floor_A": floor, "margin_above_floor_A": None if floor is None else float(minimum["distance_A"] - floor), "pass": floor is None or minimum["distance_A"] >= floor - 1e-10}
            source_ids = np.asarray(source.arrays["lammps_id"], dtype=int)
            atom_displacements = [{"id": int(source_ids[i]), "element": source[i].symbol, "delta_xyz_A": chosen["translation_A"] if source[i].symbol == "Ag" else [0.0, 0.0, 0.0]} for i in range(len(source))]
            selected_record = {
                "label": f"Ag{kind[2:]}_v16_unseen_validation_01",
                "interface": kind,
                "role": "unassigned withheld-validation geometry proposal; A must freeze and assign role",
                "status": "geometry_screen_pass_candidate_only",
                "source": {"path": str(parent_dataset_path.relative_to(repo)), "dataset_sha256": parent_sha, "frame_index_zero_based": frame_index, "source_config_type": cfg, "source_frame_geometry_sha256": meta["source_frame_geometry_sha256"], "formula": formula(source), "atoms": len(source)},
                "candidate": {"path": outpath.name, "sha256": sha(outpath), "formula": formula(check), "atoms": len(check), "cell_A": check.cell.array.tolist(), "pbc": check.pbc.tolist(), "registry_translation_A": chosen["translation_A"], "magnitude_A": chosen["magnitude_A"], "angle_from_established_axis_deg": chosen["angle_from_established_lateral_axis_deg"], "atom_displacements_by_id_A": atom_displacements, "central_pair_ids": {"Ag_id": int(pair_ids[ia]), "X_id": int(pair_ids[ix])}, "central_pair_distance_A": chosen["marked_pair_distance_A"], "observed_central_pair_distance_range_A": central_ranges[kind], "species_pair_minima": chosen["species_pair_minima"], "species_contact_screen_audit": species_contact_audit, "minimum_image_screen_pass": True, "labels_absent": True},
                "comparison": {"screened_geometry_count_same_formula": len(pools[kind]), "exact_duplicate_hits": 0, "near_duplicate_hits": 0, "nearest_known_geometry": chosen["nearest_known_geometry"], "criteria": {"exact": "same IDs/species/cell and aligned maximum displacement <=1e-5 A", "near_same_ID": f"Ag registry RMS <= {NEAR_REGISTRY_RMS_A} A and non-Ag framework RMS <= {NEAR_FRAMEWORK_RMS_A} A after ID matching and MIC", "species_contact": f"candidate species-pair minima must not compress >{PAIR_FLOOR_MARGIN_A} A below the smallest existing same-composition geometry; empirical screen, not a bond-length claim", "central_contact": "within 0.20 A of the observed same-ID marked-pair range across screened geometries"}, "source_hashes": {"historical_localization": lhash, "parent_candidate_screen_manifest": phash}, "target_framework_neighbor_changes": target_neighbor_before_after, "selection_score": score},
                "limits": ["inherits an existing small-cluster motif; not independent morphology", "does not test thermal disorder, liquid Ag or extended interface", "Ag translation changes neighboring framework environment but leaves framework coordinates fixed", "proposal has no DFT labels and is not assigned to train/validation/test until A decides"],
            }
            screening.append({"selected_candidate": selected_record, "search_grid_candidates": len(candidate_checks), "screened_grid": candidate_checks})

        # Deliver compact candidate manifest and an auditable screening record.
        all_outputs = [out / f"{k}_v16_unseen_validation_01.extxyz" for k in INTERFACE]
        candidate_manifest = {
            "title": "Window B V16 label-free unseen-validation geometry proposals",
            "owner_instance": "c5035b48-f43b-4da0-b8e4-e2862817f86a",
            "scope": {"geometry_only": True, "DFT": False, "MACE_inference_or_training": False, "MD_or_TTM": False, "V16_calculation_outputs_read": False, "V16_acquisition_inputs_read": True},
            "role_policy": "proposals are withheld-validation geometry candidates only; no split/role is assigned until A freezes and accepts them",
            "screening_policy": {"translation_plane": "xy, established z interface normal", "search_magnitudes_A": magnitudes, "angle_steps_deg": angles, "exact_duplicate_rule": "same IDs/species/cell and max aligned MIC displacement <=1e-5 A", "near_duplicate_rule": f"same IDs: Ag registry RMS <= {NEAR_REGISTRY_RMS_A} A and framework RMS <= {NEAR_FRAMEWORK_RMS_A} A; geometry matched by persistent ID and MIC", "species_contact_rule": f"for each same-composition interface, no candidate species-pair minimum may be >{PAIR_FLOOR_MARGIN_A} A below the smallest distance in the screened existing geometries; empirical guardrail, not a universal chemical bond cutoff", "central_pair_range_rule": "candidate marked-pair distance must fall within +/-0.20 A of existing observed range", "selection": "among passing candidates, maximize nearest-geometry RMS plus 0.25*maximum change in nearest-Ag distance for historically localized non-Ag target atoms"},
            "planned_v16_training_acquisition_inputs_only": planned_acquisitions,
            "sources_inventory": inventory,
            "records": [x["selected_candidate"] for x in screening],
            "input_sha256": dict(sorted(evidence_hashes.items())),
        }
        write_json(out / "candidate_manifest.json", candidate_manifest)
        write_json(out / "screening_details.json", {"selected": [x["selected_candidate"] for x in screening], "search_grid_evaluations": [{"label": x["selected_candidate"]["label"], "candidates": x["search_grid_candidates"]} for x in screening], "full_search_diagnostics": [{"selected_translation_A": x["selected_candidate"]["candidate"]["registry_translation_A"], "selected_search_score": x["selected_candidate"]["comparison"]["selection_score"], "eligible_passes": sum(1 for c in x["screened_grid"] if c["valid"]), "n_candidates": len(x["screened_grid"]), "nearest_geometry": x["selected_candidate"]["comparison"]["nearest_known_geometry"], "selected_target_changes": x["selected_candidate"]["comparison"]["target_framework_neighbor_changes"]} for x in screening]})
        report_lines = [
            "# V16 withheld-validation geometry design (label-free)",
            "",
            "## Scope",
            "",
            "Three geometries were generated from the existing V14 training parent motifs, one each for AgC, AgSi and AgTi. Only coordinates/IDs/cell/PBC, input manifests and the historical V15 force-localization report were used. No DFT calculation output, unseen DFT label, MACE model/inference/training, MD or TTM was accessed. The six planned V16 acquisition input geometries were included in deduplication; their calculation folders were not opened.",
            "",
            "The proposals move the Ag sublattice laterally while leaving the carbide/silicide framework fixed. This changes the local Ag-neighbor environment sampled by framework atoms identified in historical localization, while avoiding claims of independent morphology. Each file contains geometry/identity metadata only and remains unassigned until A freezes and assigns the split.",
            "",
            "## Selected candidates",
            "",
            "| Interface | Lateral shift (Å) | Marked pair distance (Å) | Closest known geometry RMS (Å) | Main framework-neighbor response target |",
            "|---|---:|---:|---:|---|",
        ]
        for item in candidate_manifest["records"]:
            target = item["comparison"]["target_framework_neighbor_changes"]
            target_text = ", ".join(f"{x['element']}#{x['id']} Δnearest-Ag={x['nearest_Ag_distance_change_A']:+.3f} Å" for x in target[:2]) or "none retained"
            t = item["candidate"]["registry_translation_A"]
            nearest = item["comparison"]["nearest_known_geometry"]
            report_lines.append(f"| {item['interface']} | ({t[0]:+.3f}, {t[1]:+.3f}, {t[2]:+.3f}) | {item['candidate']['central_pair_distance_A']:.3f} | {nearest['all_atom_RMS_A']:.3f} | {target_text} |")
        report_lines += [
            "",
            "The nearest-geometry, no-exact/no-near check covered every current train/validation/test frame found under all model-version data directories, prior V13–V15 frozen input geometries, the six V16 acquisition inputs, and the earlier B label-free registry/wide-span proposals. Matching uses persistent IDs where possible and same-species assignment otherwise. Full paths, hashes and nearest hits are in `candidate_manifest.json` and `screening_details.json`.",
            "",
            "## Geometry and contact audit",
            "",
            "Candidates preserve the source atom count, atomic ordering, persistent IDs, cell and PBC. The only coordinate changes are a rigid xy translation of the Ag sublattice. All pairwise element minima use minimum-image distances and are reported by species pair. The short-contact screen is comparative: candidate minima may not be more than 0.05 Å below the smallest observed distance in the screened same-composition set. It is not a universal bond-length judgment; the framework's existing C–Ti minima are preserved.",
            "",
            "No exact or near duplicate passed into the selected set. The nearest known geometry is reported for each candidate. The candidates are intentionally motif-correlated proposals, not independent-material validation. They do not test framework relaxation or thermal disorder; a separate paired framework-perturbation set would be needed to test that factor directly.",
            "",
            "## For A",
            "",
            "Treat all three as geometry proposals only. A should review the source lineage, hash/contact audit and intended validation role, then freeze accepted inputs before V16 training. Do not use them for training if their purpose remains an unseen validation set. Keep the planned six V16 acquisitions and all scored V15 probes in a distinct acquisition/history role.",
            "",
            "## Reproduction",
            "",
            "Run from the repository root in the ASE environment:",
            "",
            "```bash",
            "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\",
            "  /workspace/.venvs/mace-v14-assist/bin/python -B \\",
            "  research/mace-v12-transfer/coordination/reports/window-b/v16_unseen_validation_geometry_design/generate_candidates.py \\",
            "  --repo-root /workspace/- \\",
            "  --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v16_unseen_validation_geometry_design/reproduction_01",
            "```",
            "",
            "The script refuses nonempty output directories. The SHA inventory covers this generator, the three label-free extxyz files, manifest, screening details and this report.",
        ]
        (out / "screening_report.md").write_text("\n".join(report_lines) + "\n")
        script_copy = out / "generate_candidates.py"
        script_copy.write_bytes(Path(__file__).read_bytes())
        outputs = [script_copy, *all_outputs, out / "candidate_manifest.json", out / "screening_details.json", out / "screening_report.md"]
        (out / "SHA256SUMS.txt").write_text("".join(f"{sha(path)}  {path.name}\n" for path in outputs))
        print(json.dumps({"selected": [x["selected_candidate"] for x in screening], "grid_count": [x["search_grid_candidates"] for x in screening], "outputs": [x.name for x in outputs]}, indent=2))
    except Exception as exc:
        write_json(out / "diagnostics_failure.json", {"error": type(exc).__name__, "message": str(exc), "input_sha256": dict(sorted(evidence_hashes.items()))})
        raise


if __name__ == "__main__":
    main()
