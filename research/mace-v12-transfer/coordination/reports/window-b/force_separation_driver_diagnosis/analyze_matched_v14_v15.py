#!/usr/bin/env python3
"""Reproduce a matched V14/V15 force-driver comparison on already-scored V15 probes.

This script reads only the three archived V15 registry references already scored by
the frozen V15 evaluation. V16 references are never opened; the six V16 inputs are
read from input manifests solely for geometry/design mapping.
"""
import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path

import numpy as np

PROJECT = Path("research/mace-v12-transfer/periodic_interface_v4")
REPORT = Path("research/mace-v12-transfer/coordination/reports/window-b/force_separation_driver_diagnosis")
LABELS = [f"{x}_registry_holdout_v15_01" for x in ("AgC", "AgSi", "AgTi")]
LOCAL_CUTOFF_A = 4.5
METRIC_TOL = 1.0e-6


def need(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def metric(force_error):
    mag = np.linalg.norm(force_error, axis=1)
    return {
        "atoms": int(len(force_error)),
        "vector_RMSE_eV_A": float(np.sqrt(np.mean(mag**2))),
        "component_RMSE_xyz_eV_A": np.sqrt(np.mean(force_error**2, axis=0)).tolist(),
        "max_vector_error_eV_A": float(np.max(mag)),
        "SSE": float(np.sum(mag**2)),
    }


def region_metrics(error, symbols, mask):
    result = metric(error[mask]) if np.any(mask) else None
    if result is None:
        return None
    result.pop("SSE")
    result["species"] = {}
    for species in sorted(set(symbols[mask])):
        smask = mask & (symbols == species)
        item = metric(error[smask])
        item.pop("SSE")
        result["species"][str(species)] = item
    return result


def atoms_by_id(atoms):
    need("lammps_id" in atoms.arrays, "persistent atom IDs are missing")
    ids = np.asarray(atoms.arrays["lammps_id"], dtype=int)
    need(len(set(ids.tolist())) == len(ids), "duplicate persistent atom IDs")
    return {int(atom_id): i for i, atom_id in enumerate(ids)}


def config_type(atoms):
    return str(atoms.info.get("config_type", ""))


def mic_deltas(a, b):
    """Position differences a-b in ID order, minimum-image, before alignment."""
    from ase.geometry import find_mic

    amap, bmap = atoms_by_id(a), atoms_by_id(b)
    need(set(amap) == set(bmap), "persistent ID sets differ")
    need(np.array_equal(a.numbers, b.numbers), "atom number/order differs")
    need(np.array_equal(a.pbc, b.pbc), "PBC differs")
    need(np.allclose(a.cell.array, b.cell.array, atol=1e-9, rtol=0), "cell differs")
    d = np.asarray([a.positions[amap[k]] - b.positions[bmap[k]] for k in sorted(amap)])
    d, lengths = find_mic(d, a.cell, pbc=a.pbc)
    return sorted(amap), np.asarray(d), np.asarray(lengths)


def geometry_distance(probe, train):
    ids, delta, lengths = mic_deltas(probe, train)
    # Remove a common rigid translation using the median displacement; in these
    # isolated cluster boxes, this leaves local registry/frame changes intact.
    delta = delta - np.median(delta, axis=0)
    pmap = atoms_by_id(probe)
    ordered_symbols = np.asarray([probe[pmap[k]].symbol for k in ids])
    ag = ordered_symbols == "Ag"
    framework = ~ag
    return {
        "all_RMS_displacement_A": float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
        "Ag_RMS_displacement_A": float(np.sqrt(np.mean(np.sum(delta[ag]**2, axis=1)))) if ag.any() else None,
        "framework_RMS_displacement_A": float(np.sqrt(np.mean(np.sum(delta[framework]**2, axis=1)))) if framework.any() else None,
        "max_atom_displacement_A": float(np.max(np.linalg.norm(delta, axis=1))),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()
    repo = args.repo_root.resolve()
    out = args.output_dir.resolve()
    need(out.is_relative_to(repo / REPORT) and out != repo / REPORT, "output must be a report child")
    need(not out.exists() or (out.is_dir() and not any(out.iterdir())), "nonempty output refused")
    out.mkdir(parents=True, exist_ok=True)
    inputs = {}

    def record(path):
        path = Path(path).resolve()
        need(path.is_relative_to(repo), f"evidence escapes repository: {path}")
        need(path.is_file(), f"evidence missing: {path}")
        inputs[str(path.relative_to(repo))] = sha(path)
        return path

    def load_json(path):
        return json.loads(record(path).read_text())

    try:
        from ase.io import read
        import torch
        from mace.calculators import MACECalculator

        torch.set_default_dtype(torch.float64)
        torch.set_num_threads(1)
        version = PROJECT / "mace_periodic_v14_interface_energy"
        v15 = PROJECT / "mace_periodic_v15_interface_energy"
        eval14_path = version / "results/v14_independent_holdout_comparison.json"
        eval14_assessment_path = version / "results/v14_validation_assessment.json"
        eval15_path = v15 / "evaluation_01/evaluation.json"
        select15_path = v15 / "selected_model_record.json"
        eval14 = load_json(repo / eval14_path)
        assessment14 = load_json(repo / eval14_assessment_path)
        eval15 = load_json(repo / eval15_path)
        select15 = load_json(repo / select15_path)
        model14_path = version / "checkpoints/MACE_periodic_v14_interface_energy_run-45.model"
        model15_path = v15 / "training_01/checkpoints/MACE_periodic_v15_interface_energy_run-45.model"
        model14_file, model15_file = record(repo / model14_path), record(repo / model15_path)
        model14_sha, model15_sha = sha(model14_file), sha(model15_file)
        need(eval14["models"]["v14"]["sha256"] == model14_sha == assessment14["training"]["evaluated_model_sha256"], "V14 selected/scored model hash mismatch")
        need(select15["model_sha256"] == eval15["model_sha256"] == model15_sha, "V15 selected/scored model hash mismatch")
        need(select15["training_exit_code"] == 0 and select15["completed_epochs"] == 80 and select15["selection_uses_only_training_validation"] is True and select15["blind_labels_used_for_selection"] is False, "V15 selection evidence failed")
        need(assessment14["training"]["last_optimizer_epoch"] == 79, "V14 80-epoch completion evidence differs")
        need(assessment14["training"]["selected_checkpoint"].endswith("run-45_epoch-71.pt"), "V14 selected epoch differs from assessment")

        eval15_rows = {r["label"]: r for r in eval15["per_structure"] if r["role"] == "fresh_local_registry"}
        need(set(eval15_rows) == set(LABELS), "expected exactly the three scored V15 probes")
        manifest_path = PROJECT / "pbe_interface_v15_registry_holdouts/input_manifest.json"
        manifest = load_json(repo / manifest_path)
        meta = {r["label"]: r for r in manifest["records"]}
        need(set(meta) >= set(LABELS), "holdout input manifest is incomplete")
        train_frames = {}
        split_meta = {}
        for ver in ("v14", "v15"):
            split_meta[ver] = {}
            for split in ("train", "valid", "test"):
                path = PROJECT / f"mace_periodic_{ver}_interface_energy/data/{split}.extxyz"
                frames = read(record(repo / path), index=":")
                train_frames[(ver, split)] = frames
                split_meta[ver][split] = {"path": str(path), "sha256": inputs[str(path)], "frames": len(frames), "config_types": [config_type(a) for a in frames]}
        need(split_meta["v14"]["train"]["frames"] == 31 and split_meta["v15"]["train"]["frames"] == 29, "unexpected V14/V15 training sizes")
        for split in ("valid", "test"):
            need(split_meta["v14"][split]["sha256"] == split_meta["v15"][split]["sha256"], f"V14/V15 {split} split changed")

        # Load selected checkpoints one at a time and keep only small prediction arrays.
        predictions = {"v14": {}, "v15": {}}
        model_schema = {}
        for ver, model_file in (("v14", model14_file), ("v15", model15_file)):
            loaded = torch.load(model_file, map_location="cpu", weights_only=False).double().eval()
            model_schema[ver] = {"class": type(loaded).__name__, "state_dict_shapes": sorted((name, list(value.shape)) for name, value in loaded.state_dict().items() if torch.is_tensor(value))}
            calc = MACECalculator(models=[loaded], device="cpu", default_dtype="float64")
            for label in LABELS:
                archive = PROJECT / "pbe_interface_v15_registry_holdouts/calculations" / label
                verification = load_json(repo / archive / "verification.json")
                summary = load_json(repo / archive / "summary.json")
                progress = load_json(repo / archive / "progress.json")
                rec = meta[label]
                need(verification["status"] == "PASS" and verification["role"] == "v15_blind_registry_holdout" and all(x is True for x in verification["checks"].values()), f"archive verification failed: {label}")
                need(summary["scf_converged"] is True and progress["scf_converged"] is True and progress["status"] == "complete", f"reference not converged: {label}")
                need(verification["label"] == summary["label"] == progress["label"] == label, f"label mismatch: {label}")
                need(verification["input_sha256"] == summary["source_sha256"] == progress["source_sha256"] == rec["input_sha256"], f"input lineage mismatch: {label}")
                for member, expected_sha in verification["sha256"].items():
                    need(Path(member).name == member, "unsafe archive member")
                    need(sha(record(repo / archive / member)) == expected_sha, f"archive hash mismatch: {label}/{member}")
                refpath = archive / f"{label}_PW_PBE.extxyz"
                refpath = record(repo / refpath)
                need(sha(refpath) == eval15_rows[label]["reference_sha256"], f"reference differs from frozen evaluation: {label}")
                inpath = record(repo / (PROJECT / "pbe_interface_v15_registry_holdouts" / rec["input"]))
                need(sha(inpath) == rec["input_sha256"], f"input hash mismatch: {label}")
                probe, ref = read(inpath), read(refpath)
                need(len(probe) == len(ref) == rec["atoms"] and probe.get_chemical_formula() == ref.get_chemical_formula() == rec["formula"], f"formula/atom count mismatch: {label}")
                need(np.array_equal(probe.numbers, ref.numbers) and np.array_equal(probe.pbc, ref.pbc) and np.allclose(probe.cell.array, ref.cell.array, atol=1e-12, rtol=0) and np.allclose(probe.positions, ref.positions, atol=1e-12, rtol=0), f"geometry/cell/PBC mismatch: {label}")
                for key in ("lammps_id", "central_pair"):
                    need(key in probe.arrays and key in ref.arrays and np.array_equal(probe.arrays[key], ref.arrays[key]), f"ID/marked pair mismatch: {label}/{key}")
                need(rec["role"] == "v15_blind_registry_holdout" and rec["model_use"].startswith("Exclude from v15 training"), f"unexpected role history: {label}")
                need(label not in split_meta["v14"]["train"]["config_types"] and label not in split_meta["v15"]["train"]["config_types"], f"scored probe appears in training by label: {label}")
                ids = atoms_by_id(ref)
                need(set(ids) == set(rec["central_pair"]["persistent_ids"] + [k for k in ids if k not in rec["central_pair"]["persistent_ids"]]), "invalid ID map")
                need("PW_PBE_forces" in ref.arrays and "PW_PBE_energy_eV" in ref.info, f"reference labels missing: {label}")
                fref = np.asarray(ref.arrays["PW_PBE_forces"], dtype=float)
                eref = float(ref.info["PW_PBE_energy_eV"])
                need(fref.shape == (len(ref), 3) and np.isfinite(fref).all() and np.isfinite(eref), f"nonfinite reference: {label}")
                pred = probe.copy()
                pred.calc = calc
                energy = float(pred.get_potential_energy())
                forces = np.asarray(pred.get_forces(), dtype=float)
                need(np.isfinite(energy) and forces.shape == fref.shape and np.isfinite(forces).all(), f"nonfinite prediction: {ver}/{label}")
                predictions[ver][label] = {"energy": energy, "forces": forces, "reference_energy": eref, "reference_forces": fref, "atoms": ref, "input_sha256": rec["input_sha256"], "reference_sha256": inputs[str(refpath.relative_to(repo))]}
            del calc, loaded
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        need(model_schema["v14"] == model_schema["v15"], "V14/V15 serialized model architecture schemas differ")

        # Confirm V15 predictions reproduce the already-published scorer.
        reproduction = {}
        for label in LABELS:
            item = predictions["v15"][label]
            atoms = item["atoms"]
            symbols = np.asarray(atoms.get_chemical_symbols())
            pair = np.flatnonzero(np.asarray(atoms.arrays["central_pair"], dtype=bool))
            need(len(pair) == 2, f"marked pair missing: {label}")
            ia = next(int(i) for i in pair if symbols[i] == "Ag")
            ix = next(int(i) for i in pair if symbols[i] != "Ag")
            vec = atoms.positions[ix] - atoms.positions[ia]
            unit = vec / np.linalg.norm(vec)
            radial_ref = float(np.dot(item["reference_forces"][ix] - item["reference_forces"][ia], unit))
            radial_pred = float(np.dot(item["forces"][ix] - item["forces"][ia], unit))
            ferr = item["forces"] - item["reference_forces"]
            reproduced = {
                "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(np.sum(ferr**2, axis=1)))),
                "force_vector_max_error_eV_A": float(np.linalg.norm(ferr, axis=1).max()),
                "separating_force_abs_error_eV_A": abs(radial_pred - radial_ref),
                "energy_prediction_eV_cell": item["energy"],
            }
            target = eval15_rows[label]
            diffs = {k: abs(reproduced[k] - target[k]) for k in reproduced}
            need(all(v <= METRIC_TOL for v in diffs.values()), f"V15 published metric reproduction failed: {label}: {diffs}")
            reproduction[label] = {"reproduced": reproduced, "published": {k: target[k] for k in reproduced}, "abs_deltas": diffs}

        rows = []
        decompositions = {"input_scope": "three already-scored V15 probes only; no V16 reference/output opened", "local_cutoff_A": LOCAL_CUTOFF_A, "structures": []}
        parent_similarity = {}
        for label in LABELS:
            ref_item = predictions["v14"][label]
            atoms = ref_item["atoms"]
            symbols = np.asarray(atoms.get_chemical_symbols())
            ids = np.asarray(atoms.arrays["lammps_id"], dtype=int)
            dist = atoms.get_all_distances(mic=True)
            pair = np.flatnonzero(np.asarray(atoms.arrays["central_pair"], dtype=bool))
            need(len(pair) == 2, f"marked pair missing: {label}")
            ia = next(int(i) for i in pair if symbols[i] == "Ag")
            ix = next(int(i) for i in pair if symbols[i] != "Ag")
            need(set(ids[pair].tolist()) == set(meta[label]["central_pair"]["persistent_ids"]), f"central IDs differ from manifest: {label}")
            direct = atoms.positions[ix] - atoms.positions[ia]
            shortest = atoms.get_distance(ia, ix, mic=True, vector=True)
            direct_norm = float(np.linalg.norm(direct))
            mic_norm = float(np.linalg.norm(shortest))
            need(abs(direct_norm - mic_norm) < 1e-8, f"direct/MIC pair image ambiguity: {label}")
            u = direct / direct_norm
            central_local_distance = np.min(dist[:, pair], axis=1)
            local = central_local_distance <= LOCAL_CUTOFF_A
            marked = np.zeros(len(atoms), dtype=bool); marked[pair] = True
            ag_neighbors = local & (symbols == "Ag") & ~marked
            framework_shell = local & (symbols != "Ag") & ~marked
            outside = ~local
            categories = {"marked_pair": marked, "Ag_neighbors_local": ag_neighbors, "framework_shell_local": framework_shell, "outside_4p5A": outside}
            need(sum(int(x.sum()) for x in categories.values()) == len(atoms), "region masks do not partition all atoms")
            model_data = {}
            for ver in ("v14", "v15"):
                entry = predictions[ver][label]
                error = entry["forces"] - entry["reference_forces"]
                radial_ref = float(np.dot(entry["reference_forces"][ix] - entry["reference_forces"][ia], u))
                radial_pred = float(np.dot(entry["forces"][ix] - entry["forces"][ia], u))
                contributions = []
                for index, sign in ((ia, -1.0), (ix, 1.0)):
                    projection = float(np.dot(error[index], u))
                    trans = error[index] - projection * u
                    contributions.append({
                        "atom_id": int(ids[index]), "element": str(symbols[index]),
                        "force_error_xyz_eV_A": error[index].tolist(),
                        "longitudinal_error_projection_eV_A": projection,
                        "signed_contribution_to_separation_error_eV_A": sign * projection,
                        "transverse_error_xyz_eV_A": trans.tolist(),
                        "transverse_error_norm_eV_A": float(np.linalg.norm(trans)),
                        "full_vector_error_norm_eV_A": float(np.linalg.norm(error[index])),
                    })
                separation_error = radial_pred - radial_ref
                need(abs(sum(x["signed_contribution_to_separation_error_eV_A"] for x in contributions) - separation_error) < 1e-10, "separation projection terms do not sum")
                mag = np.linalg.norm(error, axis=1)
                top = []
                for i in np.argsort(mag)[::-1][:10]:
                    js = [j for j in np.argsort(dist[i]) if j != i][:4]
                    top.append({"atom_id": int(ids[i]), "element": str(symbols[i]), "force_error_norm_eV_A": float(mag[i]), "force_error_xyz_eV_A": error[i].tolist(), "distance_to_marked_pair_min_MIC_A": float(central_local_distance[i]), "region": next(k for k, m in categories.items() if m[i]), "nearest4_MIC": [{"atom_id": int(ids[j]), "element": str(symbols[j]), "distance_A": float(dist[i, j])} for j in js]})
                regions = {name: region_metrics(error, symbols, mask) for name, mask in categories.items()}
                for name, mask in categories.items():
                    if regions[name] is not None:
                        regions[name]["force_error_SSE_share"] = float(np.sum(mag[mask]**2) / np.sum(mag**2)) if np.sum(mag**2) else 0.0
                model_data[ver] = {
                    "model_sha256": model14_sha if ver == "v14" else model15_sha,
                    "energy_prediction_eV_cell": entry["energy"],
                    "energy_error_meV_atom_signed": float((entry["energy"] - entry["reference_energy"]) * 1000 / len(atoms)),
                    "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(mag**2))),
                    "force_component_RMSE_xyz_eV_A": np.sqrt(np.mean(error**2, axis=0)).tolist(),
                    "force_vector_max_error_eV_A": float(mag.max()),
                    "separation_reference_eV_A": radial_ref,
                    "separation_prediction_eV_A": radial_pred,
                    "separation_error_signed_eV_A": separation_error,
                    "separation_error_abs_eV_A": abs(separation_error),
                    "marked_pair_error_projection": contributions,
                    "error_regions": regions,
                    "top10_atom_errors": top,
                }
                row = {"label": label, "interface": label.split("_")[0], "model": ver, "atoms": len(atoms), "formula": atoms.get_chemical_formula(), "reference_sha256": entry["reference_sha256"], **{k: model_data[ver][k] for k in ("energy_error_meV_atom_signed", "force_vector_RMSE_eV_A", "force_component_RMSE_xyz_eV_A", "force_vector_max_error_eV_A", "separation_reference_eV_A", "separation_prediction_eV_A", "separation_error_signed_eV_A", "separation_error_abs_eV_A")}}
                rows.append(row)
            decompositions["structures"].append({
                "label": label, "formula": atoms.get_chemical_formula(), "atoms": len(atoms),
                "marked_pair": {"Ag_id": int(ids[ia]), "X_id": int(ids[ix]), "X_element": str(symbols[ix]), "direct_distance_A": direct_norm, "MIC_distance_A": mic_norm, "unit_Ag_to_X_direct": u.tolist()},
                "roles": {"evaluation_role_at_scoring": "V15 fresh_local_registry scored post-freeze; acquisition/history only now", "source_config_type": meta[label]["source_config_type"], "source_dataset_sha256": meta[label]["source_dataset_sha256"], "source_frame_geometry_sha256": meta[label]["source_frame_geometry_sha256"], "parent_motif_note": meta[label]["independence_limit"]},
                "local_shell_rule": "minimum-image distance <=4.5 A from either marked atom; same historical localization cutoff, not fitted to these scores",
                "models": model_data,
            })

            # Compare each probe against same persistent-ID training parents; structural only.
            parent = meta[label]["source_config_type"]
            parent_similarity[label] = {}
            for ver in ("v14", "v15"):
                candidates = []
                for train in train_frames[(ver, "train")]:
                    if config_type(train) != parent:
                        continue
                    d = geometry_distance(atoms, train)
                    candidates.append({"config_type": config_type(train), **d})
                need(candidates, f"source parent absent from {ver} train: {label}")
                parent_similarity[label][ver] = {"parent_in_train": True, "source_parent_comparison": candidates[0]}

        # Compare split membership and exact dataset edits.
        cfg14 = set(split_meta["v14"]["train"]["config_types"])
        cfg15 = set(split_meta["v15"]["train"]["config_types"])
        removed, added = sorted(cfg14 - cfg15), sorted(cfg15 - cfg14)
        need(len(removed) == 5 and len(added) == 3, f"unexpected V14->V15 train edits: removed={removed}, added={added}")
        need(all("LCAO" in x or "2p387" in x for x in removed), "unexpected excluded V14 train record")
        need(all("residual_shell" in x for x in added), "unexpected V15 additions")
        removed_rows = []
        for cfg in removed:
            frame = next(a for a in train_frames[("v14", "train")] if config_type(a) == cfg)
            removed_rows.append({"config_type": cfg, "source_method": frame.info.get("source_method"), "dft_method": frame.info.get("dft_method"), "energy_present": "REF_energy" in frame.info, "forces_present": "REF_forces" in frame.arrays, "evidence_status": "four force-only LCAO frames removed" if "LCAO" in cfg else "AgTi original source unresolved; excluded"})
        added_rows = []
        for cfg in added:
            frame = next(a for a in train_frames[("v15", "train")] if config_type(a) == cfg)
            added_rows.append({"config_type": cfg, "source_method": frame.info.get("source_method"), "dft_method": frame.info.get("dft_method"), "energy_present": "REF_energy" in frame.info, "forces_present": "REF_forces" in frame.arrays, "source_label": cfg.split("_periodic")[0]})

        # Validate and characterize the six planned V16 geometries without touching calculations/outputs.
        planned_v16 = []
        for dirname in ("pbe_interface_v16_parallel_acquisition", "pbe_interface_v16_main_acquisition"):
            manifest_path = PROJECT / dirname / "input_manifest.json"
            d = load_json(repo / manifest_path)
            need(len(d["records"]) == 3, f"expected three planned geometries in {dirname}")
            for rec in d["records"]:
                ipath = PROJECT / dirname / rec["input"]
                fp = record(repo / ipath)
                need(sha(fp) == rec["input_sha256"], f"V16 planned input hash mismatch: {rec['label']}")
                geom = read(fp)
                need(len(geom) == rec["atoms"] and geom.get_chemical_formula() == rec["formula"], f"V16 geometry identity mismatch: {rec['label']}")
                need("REF_energy" not in geom.info and "PW_PBE_energy_eV" not in geom.info and not any(k in geom.arrays for k in ("REF_forces", "PW_PBE_forces")), f"V16 input contains labels: {rec['label']}")
                central = np.flatnonzero(np.asarray(geom.arrays["central_pair"], dtype=bool))
                need(len(central) == 2 and set(np.asarray(geom.arrays["lammps_id"], dtype=int)[central].tolist()) == set(rec["central_pair"]["persistent_ids"]), f"V16 marked IDs mismatch: {rec['label']}")
                pair_distance = float(geom.get_distance(*central, mic=False))
                parent = next((a for a in train_frames[("v14", "train")] if config_type(a) == rec["parent_lineage"]["config_type"]), None)
                need(parent is not None, f"planned V16 parent not found in V14 train: {rec['label']}")
                pmap, gmap = atoms_by_id(parent), atoms_by_id(geom)
                need(set(pmap) == set(gmap), f"planned V16 parent ID set differs: {rec['label']}")
                pair_ids = rec["central_pair"]["persistent_ids"]
                ag_id = int(rec["central_pair"]["Ag_id"])
                x_id = int(rec["central_pair"]["contact_id"])
                parent_vec = parent.positions[pmap[x_id]] - parent.positions[pmap[ag_id]]
                candidate_vec = geom.positions[gmap[x_id]] - geom.positions[gmap[ag_id]]
                pair_vector_delta = candidate_vec - parent_vec
                parent_pair_distance = float(np.linalg.norm(parent_vec))
                need(abs(parent_pair_distance - float(rec["parent_lineage"].get("distance_A", parent_pair_distance))) < 1e-6, f"planned V16 parent pair distance changed from manifest: {rec['label']}")
                # The parent-frame vector delta distinguishes normal distance changes from a registry shift.
                planned_v16.append({"label": rec["label"], "candidate_id": rec.get("candidate_id"), "owner": rec.get("proposed_owner"), "role": rec["role"], "input_sha256": rec["input_sha256"], "atoms": len(geom), "formula": geom.get_chemical_formula(), "central_pair": rec["central_pair"], "parent_direct_distance_A": parent_pair_distance, "candidate_direct_distance_A": pair_distance, "direct_distance_delta_A": float(pair_distance - parent_pair_distance), "Ag_to_X_pair_vector_delta_from_parent_A": pair_vector_delta.tolist(), "parent_config_type": rec["parent_lineage"]["config_type"], "parent_geometry_sha256": rec["parent_lineage"]["geometry_sha256"], "design_axis": "radial/distance" if "parallel" in dirname else "lateral registry", "can_test": "Ag-X normal distance response" if "parallel" in dirname else "lateral Ag registry/contact response", "cannot_test": "independent framework morphology/thermal response; provenance of inherited labels; optimizer causality"})

        # Training/evaluation settings and checkpoint-selection evidence.
        run14 = record(repo / version / "run_v14_training.sh")
        run15 = record(repo / v15 / "run_v15_training.sh")
        foundation_path = record(repo / PROJECT / "mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model")
        foundation_sha = sha(foundation_path)
        run14_text, run15_text = run14.read_text(), run15.read_text()
        common_cli = {"seed": "45", "device": "cpu", "batch_size": "1", "max_num_epochs": "80", "lr": "0.0001", "weight_decay": "5e-7", "energy_weight": "100", "forces_weight": "1000"}
        for flag, value in common_cli.items():
            need(f"--{flag} {value}" in run14_text and f"--{flag} {value}" in run15_text, f"V14/V15 declared training option differs: --{flag}")
        trainlog14 = record(repo / version / "results/MACE_periodic_v14_interface_energy_run-45_train.txt")
        trainlog15 = record(repo / v15 / "training_01/results/MACE_periodic_v15_interface_energy_run-45_train.txt")
        def parse_train_log(path):
            parsed = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
            evals = [x for x in parsed if x.get("mode") == "eval" and x.get("epoch") is not None]
            need(len(evals) == 80 and {int(x["epoch"]) for x in evals} == set(range(80)), f"training eval curve incomplete: {path}")
            return evals
        log14, log15 = parse_train_log(trainlog14), parse_train_log(trainlog15)
        selected14 = 70  # one-based checkpoint epoch 71
        selected15 = int(select15["selected_epoch_zero_based"])
        training_metrics = {}
        for ver, logs, selected in (("v14", log14, selected14), ("v15", log15, selected15)):
            need(0 <= selected < 80, f"selected epoch invalid: {ver}")
            best = min(logs, key=lambda x: x["rmse_f"])
            training_metrics[ver] = {"epochs_evaluated": len(logs), "selected_epoch_zero_based": selected, "selected_validation": {k: logs[selected][k] for k in ("mae_e_per_atom", "rmse_e_per_atom", "mae_f", "rmse_f")}, "minimum_validation_force_RMSE": {"epoch_zero_based": int(best["epoch"]), "rmse_f": best["rmse_f"]}, "epoch79_validation": {k: logs[-1][k] for k in ("mae_e_per_atom", "rmse_e_per_atom", "mae_f", "rmse_f")}}

        # Build common-score CSV with explicit V15-V14 differences.
        csv_by_label = {label: {r["model"]: r for r in rows if r["label"] == label} for label in LABELS}
        common_rows = []
        for label in LABELS:
            r14, r15 = csv_by_label[label]["v14"], csv_by_label[label]["v15"]
            common_rows.append({"label": label, "interface": r14["interface"], "atoms": r14["atoms"], "formula": r14["formula"], "reference_sha256": r14["reference_sha256"],
                "v14_energy_abs_error_meV_atom": abs(r14["energy_error_meV_atom_signed"]), "v15_energy_abs_error_meV_atom": abs(r15["energy_error_meV_atom_signed"]), "delta_v15_minus_v14_energy_abs_error_meV_atom": abs(r15["energy_error_meV_atom_signed"]) - abs(r14["energy_error_meV_atom_signed"]),
                "v14_force_vector_RMSE_eV_A": r14["force_vector_RMSE_eV_A"], "v15_force_vector_RMSE_eV_A": r15["force_vector_RMSE_eV_A"], "delta_v15_minus_v14_force_vector_RMSE_eV_A": r15["force_vector_RMSE_eV_A"] - r14["force_vector_RMSE_eV_A"],
                "v14_force_component_RMSE_xyz_eV_A": json.dumps(r14["force_component_RMSE_xyz_eV_A"], separators=(",", ":")), "v15_force_component_RMSE_xyz_eV_A": json.dumps(r15["force_component_RMSE_xyz_eV_A"], separators=(",", ":")),
                "v14_separation_abs_error_eV_A": r14["separation_error_abs_eV_A"], "v15_separation_abs_error_eV_A": r15["separation_error_abs_eV_A"], "delta_v15_minus_v14_separation_abs_error_eV_A": r15["separation_error_abs_eV_A"] - r14["separation_error_abs_eV_A"],
                "v14_max_vector_error_eV_A": r14["force_vector_max_error_eV_A"], "v15_max_vector_error_eV_A": r15["force_vector_max_error_eV_A"]})
        with (out / "common_geometry_scores.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(common_rows[0]), lineterminator="\n")
            w.writeheader(); w.writerows(common_rows)
        write_json(out / "force_decomposition.json", decompositions)

        # Original source declarations / conventions are correlated with dataset inclusion, not confirmed as causes.
        provenance_path = repo / "research/mace-v12-transfer/coordination/reports/window-b/v15_inherited_method_provenance_trace/per_frame_ledger.json"
        provenance = load_json(provenance_path)
        method_review_path = repo / "research/mace-v12-transfer/coordination/reports/window-b/v15_reference_convention_review/convention_review.json"
        method_review = load_json(method_review_path)
        evidence = {
            "models": {"v14": {"sha256": model14_sha, "path": str(model14_path), "selected_checkpoint_one_based_epoch": 71, "evidence_path": str(eval14_assessment_path), "status": assessment14["status"]}, "v15": {"sha256": model15_sha, "path": str(model15_path), "selection_sha256": inputs[str(select15_path)], "selected_checkpoint_one_based_epoch": select15["selected_epoch"], "selection_uses_only_training_validation": select15["selection_uses_only_training_validation"], "status": eval15["overall_status"]}},
            "reference_and_geometry": {"probe_labels": LABELS, "all_probes_already_scored_in_v15": True, "v15_published_score_reproduction": reproduction, "localization_cutoff_A": LOCAL_CUTOFF_A, "direct_cluster_vs_MIC_pair_image_checked": True, "reference_order": "output and frozen input have exact same numbers, lammps_id, central_pair, positions, cell, PBC; reference force vectors are read at matching serialized atom indices"},
            "splits": split_meta,
            "v14_to_v15_train_edit": {"v14_frames": 31, "v15_frames": 29, "removed": removed_rows, "added": added_rows, "validation_and_test_byte_identical": True, "holdout_input_labels_excluded_from_both_training_splits": True, "history": "V15 probe geometries are translated registry variants of existing V14 training motifs; they were post-freeze V15 evaluation and are now scored acquisition/history, not blind data."},
            "inherited_provenance": {"source_path": str(provenance_path.relative_to(repo)), "reported_counts": provenance["counts"], "fully_confirmed_original_calculations": provenance["fully_confirmed_original_calculations"], "reviewed_method_convention": {"source_path": str(method_review_path.relative_to(repo)), "summary": "14/15 historical recovery rows have partial method declarations without original run/convergence confirmation; 1 unresolved AgTi row was excluded. Physical method consistency remains UNKNOWN. Three V15 residual labels have stronger direct runner/log/archive linkage; no inference that their method is historically universal."}},
            "training": {"cli_settings_confirmed_by_hashed_run_scripts": {"foundation_model": "mace-mp-0b3-medium.model", "foundation_model_sha256": foundation_sha, **common_cli, "energy_key": "REF_energy", "forces_key": "REF_forces", "selected_by": "training/validation only"}, "model_architecture_state_dict_schema_identical": True, "model_architecture_schema": model_schema["v14"], "validation_curves": training_metrics, "v14_run_script_sha256": inputs[str(run14.relative_to(repo))], "v15_run_script_sha256": inputs[str(run15.relative_to(repo))]},
            "reference_driver": {"runner_path": "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v15_registry_holdouts/run_pw_registry_holdout.py", "runner_sha256": inputs[str(record(repo / "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v15_registry_holdouts/run_pw_registry_holdout.py").relative_to(repo))], "summary_method": "fixed geometry GPAW PW-PBE 500 eV, Gamma 1x1x1, FermiDirac 0.1 eV; energy is native extrapolated from get_potential_energy(force_consistent=False), free energy separately stored; get_forces separately requested; SCF thresholds 1e-5; archived force array ordered with same ASE atoms object and checked against frozen input", "known_limit": "same run recipe is directly documented for these three labels only; historical inherited frames do not all have original logs/getter/convergence evidence"},
            "rmse_definition": {"vector": "sqrt(mean_i(sum_xyz((F_pred-F_ref)^2)))", "component": "sqrt(mean_all_atoms_and_xyz((F_pred-F_ref)^2))", "ratio": "vector RMSE = sqrt(3) times component RMSE for a single complete atom set", "separation": "(F_X-F_Ag) dot direct cluster-coordinate unit vector; signed Ag and X residual projections are reported separately", "PBC": "MIC for nearest-neighbor/local-shell distances; legacy separation score uses direct cluster vector, and direct vs MIC pair distance is verified equal for these probes"},
            "v16_inputs_only": planned_v16,
            "input_sha256": dict(sorted(inputs.items())),
            "runtime": {"python": os.sys.version.split()[0], "torch": torch.__version__, "mace_torch": importlib.metadata.version("mace-torch"), "ase": importlib.metadata.version("ase"), "device": "CPU", "torch_threads": torch.get_num_threads()},
        }
        write_json(out / "audit_evidence.json", evidence)

        # Concise ranked evidence table is emitted as a standalone reviewable report.
        report = f"""# V14/V15 force and separating-force driver diagnosis

## Scope and status

This is a matched retrospective comparison: both frozen models were evaluated on the same three already-scored V15 registry probes. CPU float64 inference used V14 model `{model14_sha}` and V15 model `{model15_sha}`. The exact V15 scores reproduce the frozen evaluation within {METRIC_TOL:g} eV/Å or eV/cell. All probe input/reference hashes, IDs, formula, coordinates, cell, PBC, and SCF verification were checked. The three V15 geometries are no longer blind; their labels are acquisition/history. No V16 result/reference, new DFT, training, MD, or TTM was used.

## Common geometry comparison

| Interface | V14 force vector RMSE | V15 force vector RMSE | V15−V14 | V14 separation error | V15 separation error | V15−V14 |
|---|---:|---:|---:|---:|---:|---:|
""" + "\n".join(f"| {r['interface']} | {r['v14_force_vector_RMSE_eV_A']:.6f} | {r['v15_force_vector_RMSE_eV_A']:.6f} | {r['delta_v15_minus_v14_force_vector_RMSE_eV_A']:+.6f} | {r['v14_separation_abs_error_eV_A']:.6f} | {r['v15_separation_abs_error_eV_A']:.6f} | {r['delta_v15_minus_v14_separation_abs_error_eV_A']:+.6f} |" for r in common_rows) + f"""

Energy MAE remains small on these geometries; exact model-by-model energy, component RMSE and max-vector scores are in `common_geometry_scores.csv`. Deltas are descriptive on a shared input, not causal estimates. The probes share small-cluster source motifs that appear in both training sets; this is a useful controlled geometry comparison but weak evidence for transfer to a different morphology.

## What changes on the same structures

""" + "\n".join(f"- **{r['interface']}:** force RMSE {r['v14_force_vector_RMSE_eV_A']:.4f} → {r['v15_force_vector_RMSE_eV_A']:.4f} eV/Å; separation absolute error {r['v14_separation_abs_error_eV_A']:.4f} → {r['v15_separation_abs_error_eV_A']:.4f} eV/Å." for r in common_rows) + f"""

V14→V15 train edits were five removals and three additions: four force-only LCAO rows and one unresolved AgTi row left; three directly archived residual-shell labels entered. V15 train size is 29 versus V14's 31. Validation and test files are byte-identical. V15 preserved inherited method declarations where available, but the provenance review found zero fully confirmed inherited original calculations among its 15 recovery targets; 14 are partial and one unresolved row was excluded. V15's residual references have direct runner/log/archive evidence. These dataset/provenance edits happened together, so the matched delta cannot identify which edit helped or hurt.

Both training runs used the same declared CLI settings: seed 45, CPU, batch size 1, 80 epochs, learning rate 1e-4, weight decay 5e-7, energy weight 100 and force weight 1000, with the same foundation-model SHA. The selected models also have matching serialized state-dictionary key/shape schemas. Checkpoints were selected from training/validation, V14 one-based epoch 71 and V15 one-based epoch {select15['selected_epoch']}. V15's selected validation force RMSE was {training_metrics['v15']['selected_validation']['rmse_f']:.5f} eV/Å; the V14 selected value was {training_metrics['v14']['selected_validation']['rmse_f']:.5f}. At epoch 79 those values were {training_metrics['v14']['epoch79_validation']['rmse_f']:.5f} and {training_metrics['v15']['epoch79_validation']['rmse_f']:.5f}. Validation improved slightly, while this matched local registry set still fails force gates; validation does not represent every local transverse response.

## Spatial and projection diagnosis

Forces are partitioned without overlap into the marked Ag/X pair, other Ag atoms within 4.5 Å MIC of either marked atom, non-Ag framework atoms in that same shell, and all atoms outside it. This 4.5 Å cutoff matches the previous localization analysis; it was not selected by optimizing these results. Each category includes vector and x/y/z component RMSE, species breakdown and SSE share in `force_decomposition.json`; top-error atoms include persistent IDs and nearest MIC neighbors.

        """

        # The markdown above is augmented below using the computed full decompositions.
        region_lines = []
        for s in decompositions["structures"]:
            v14d, v15d = s["models"]["v14"], s["models"]["v15"]
            pair14 = v14d["marked_pair_error_projection"]
            region_lines.append(f"### {s['label'].split('_')[0]} (Ag#{s['marked_pair']['Ag_id']}–{s['marked_pair']['X_element']}#{s['marked_pair']['X_id']})")
            region_lines.append("")
            region_lines.append("| Model | Pair Ag signed term | Pair X signed term | Sum separation error | Ag transverse norm | X transverse norm | Local Ag-neighbor RMSE | Local framework RMSE | Outside RMSE |")
            region_lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
            for ver in ("v14", "v15"):
                m = s["models"][ver]; terms = m["marked_pair_error_projection"]
                agterm = next(t for t in terms if t["element"] == "Ag")
                xterm = next(t for t in terms if t["element"] == s["marked_pair"]["X_element"])
                reg = m["error_regions"]
                ag_rmse = reg["Ag_neighbors_local"]["vector_RMSE_eV_A"] if reg["Ag_neighbors_local"] else None
                fw_rmse = reg["framework_shell_local"]["vector_RMSE_eV_A"] if reg["framework_shell_local"] else None
                out_rmse = reg["outside_4p5A"]["vector_RMSE_eV_A"] if reg["outside_4p5A"] else None
                region_lines.append(f"| {ver} | {agterm['signed_contribution_to_separation_error_eV_A']:+.4f} | {xterm['signed_contribution_to_separation_error_eV_A']:+.4f} | {m['separation_error_signed_eV_A']:+.4f} | {agterm['transverse_error_norm_eV_A']:.4f} | {xterm['transverse_error_norm_eV_A']:.4f} | {ag_rmse if ag_rmse is not None else float('nan'):.4f} | {fw_rmse if fw_rmse is not None else float('nan'):.4f} | {out_rmse if out_rmse is not None else float('nan'):.4f} |")
            region_lines.append("")
            highest = s["models"]["v15"]["top10_atom_errors"][:5]
            region_lines.append("V15 top five atom errors: " + "; ".join(f"{a['element']}#{a['atom_id']} {a['force_error_norm_eV_A']:.3f} eV/Å" for a in highest) + ".")
            region_lines.append("")
        report = report.rstrip() + "\n\n" + "\n".join(region_lines)
        report += """

The marked Ag and X terms sum exactly to the scalar separation error. Their transverse residuals are not included in that scalar. A small radial error can therefore coexist with a large vector error; changing only distance or radial-force weighting would not target that failure. Conversely, a large scalar separation error can have both pair terms aligned in the same direction. Pair distance was checked under both direct cluster coordinates and MIC; those agree for all three probes. Local neighbor distances use MIC, while the legacy signed separation score uses the direct Ag→X vector.

## Ranked hypotheses and decisions for A

| Rank | Evidence and counterevidence | Confidence | Confound | Discriminating test | Suggested action |
|---|---|---|---|---|---|
| 1. Registry/transverse coverage and local environment | Same-geometry decomposition identifies local Ag/framework errors, including transverse components; AgSi's scalar separation can look small while its vector error remains large. The three probes share inherited cluster motifs and only sample a few registry shifts. | Medium-high that this is a concrete error mode; low on its share of the V15→V14 change. | V15 simultaneously changed labels and exclusions; only one geometry per interface here. | On fixed parent structures, compare paired ± in-plane shifts at constant Ag–X normal distance; separately perturb framework coordinates with Ag held fixed. Keep identical DFT method and training recipe. | V16's three registry inputs probe lateral registry, while its three distance inputs probe normal response. They do not test independent framework motion or new morphology. A should review all six labels before deciding inclusion. |
| 2. Historical reference-method uncertainty | Provenance review documents 14 partially recovered inherited labels and one unresolved AgTi row excluded; inherited energy conventions and convergence are not uniformly proven. The three current probes themselves have passing verification and direct matched GPAW archives. | Medium that this is a dataset-risk factor; low that it explains these specific force residuals. | The evaluated probes were newly calculated with the same documented protocol; no paired re-label data exists. | Recover source inputs/logs, or make a separately authorized exact-geometry matched-method reference set with full getter/energy convention and force-order metadata. | Keep current labels unchanged; do not retrofit free energy or infer missing setup. Preserve per-row evidence tiers and isolate convention comparisons. |
| 3. Training objective/optimization or checkpoint generalization | V15 validation force RMSE is modestly better than V14 but held-out registry force errors do not meet the gate. Both runs used same declared settings; no evidence singles out loss weights, learning rate, epoch count or architecture. | Low to medium. | Changed data composition, inherited target quality and validation coverage; single seed and tiny validation split. | After a fixed data/reference baseline is frozen, run one-factor matched ablations (first force-weight 1000 vs 2000 with seed/settings held fixed; select only on a strengthened validation set). | For A's immediate V16 cycle, preserve the current CLI settings so the six new structures' data effect can be read. Do not increase epochs or force weight as an untested repair. |
| 4. Force/vector scoring, atom order or PBC bug | The frozen V15 scores are reproduced; model/reference structures match atom IDs, order, positions, cell and PBC. Vector RMSE and Cartesian component RMSE are separately reported. Pair separation uses the documented direct vector, which agrees with MIC distance for these probes. | Low for these three comparisons; checks rule out common implementation/order errors. | This audit validates these exact archives and scorer convention, not every inherited dataset serialization or every future structure. | Keep automated exact-ID/order and vector-vs-component regression checks in future evaluation. | Retain the existing metric definition; do not replace it with a component RMSE or change the PBC convention silently. |
| 5. DFT energy/smearing/k-point mismatch as direct force cause | Current probes use a documented PBE plane-wave/Gamma/Fermi-Dirac 0.1 eV driver; native extrapolated energy and free energy are separately stored. Historical labels show multiple cells/sampling conventions and incomplete original evidence. Energy errors on current probes are small. | Low for a direct cause on these references; broader target consistency remains unresolved. | Energy offsets do not imply force errors, and no matched convention perturbation exists. | Only if provenance recovery warrants it, compare native/free energy and forces on exact same geometries under documented converged conditions. | Keep REF_energy convention fixed and explicit; record both energies, smearing, k-points, forces, software and convergence for future labels. |

The parameter recommendation is intentionally conservative and directly usable: for the next data-integration cycle, retain seed 45, CPU, batch 1, 80 epochs, lr 1e-4, weight decay 5e-7, energy weight 100 and force weight 1000 while adding only reviewed V16 labels. That keeps the training recipe controlled so A can attribute differences more cleanly to the changed training set. If errors remain after source and split audits, a future fixed-data two-arm force-weight comparison (1000 vs 2000) is a proposal, not an instruction to train now; use additional validation geometries for checkpoint selection and keep a new independent blind set untouched.

## Six planned V16 inputs: what they test

""" + "\n".join(f"- `{x['label']}` ({x['candidate_id']}, {x['design_axis']}, direct gap {x['candidate_direct_distance_A']:.3f} Å; parent {x['parent_direct_distance_A']:.3f} Å): {x['can_test']}; it cannot resolve {x['cannot_test']}. Input hash `{x['input_sha256']}`." for x in planned_v16) + """

The six are training acquisitions, not blind validation. The parallel-set distances are reported in `audit_evidence.json` against the parent pair value; the main-set structures primarily alter registry. The minimal follow-up proposal is paired ±0.15 Å transverse translations (two independent in-plane axes, fixed normal separation) and paired ±0.04 Å local framework displacements for selected C/Ti/Si neighbors, with the Ag sublattice fixed for the framework arm. Generate matching controls from the same parent and use one unchanged PBE/PW/smearing/k-point/convergence recipe. Assign labels only after A approves. Keep any new independent-morphology/thermal validation inputs separate from these acquisition variants.

## Reproducibility and limits

Run from the repository root with the B MACE CPU environment:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
  /workspace/.venvs/mace-v14-assist/bin/python -B \\
  research/mace-v12-transfer/coordination/reports/window-b/force_separation_driver_diagnosis/analyze_matched_v14_v15.py \\
  --repo-root /workspace/- \\
  --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/force_separation_driver_diagnosis/analysis_rerun_01
```

Outputs: `common_geometry_scores.csv`, `force_decomposition.json`, `audit_evidence.json`, this report and `SHA256SUMS.txt`. Evidence hashes cover both selected models, selected/evaluation records, three scored references and verifications, source input manifests/frames, V14/V15 split and training logs/scripts, inherited provenance/convention audits, V16 input manifests/frames, and the exact runner. V16 calculation directories are not read. No new data labels, model fitting, or model selection were produced.
"""
        (out / "diagnosis_report.md").write_text(report)
        write_json(out / "v15_score_reproduction.json", reproduction)
        # Evidence file now includes the hashes acquired by record() through the full analysis.
        evidence["input_sha256"] = dict(sorted(inputs.items()))
        write_json(out / "audit_evidence.json", evidence)
        script_copy = out / "analyze_matched_v14_v15.py"
        script_copy.write_bytes(Path(__file__).read_bytes())
        inventory_paths = [out / p for p in ("analyze_matched_v14_v15.py", "common_geometry_scores.csv", "force_decomposition.json", "audit_evidence.json", "v15_score_reproduction.json", "diagnosis_report.md")]
        (out / "SHA256SUMS.txt").write_text("".join(f"{sha(p)}  {p.name}\n" for p in inventory_paths))
        print(json.dumps({"common_geometry_scores": common_rows, "v15_score_reproduction": reproduction, "outputs": [p.name for p in inventory_paths]}, indent=2))
    except Exception as exc:
        write_json(out / "diagnostics_failure.json", {"error": type(exc).__name__, "message": str(exc), "input_sha256": dict(sorted(inputs.items()))})
        raise


if __name__ == "__main__":
    main()
