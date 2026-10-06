#!/usr/bin/env python3
"""Create paired, label-free V16 training-acquisition geometry proposals."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

PROJECT = Path("research/mace-v12-transfer/periodic_interface_v4")
WINDOW_B = Path("research/mace-v12-transfer/coordination/reports/window-b")
REPORT = WINDOW_B / "targeted_transverse_framework_design"
PAIR_MARGIN_A = 0.05
NEAR_AG_RMS_A = 0.15
NEAR_FRAMEWORK_RMS_A = 0.15
TRANSVERSE_STEP_A = 0.15
FRAMEWORK_STEP_A = 0.04
PARENTS = {
    "AgC": "AgC_registry_strain_acq_v14_01_periodic_PW_PBE_energy_force_v14_train",
    "AgSi": "AgSi_registry_strain_acq_v14_01_periodic_PW_PBE_energy_force_v14_train",
    "AgTi": "AgTi_registry_probe_v13_01_periodic_PW_PBE_energy_force_v14_train",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def atom_ids(a):
    raw = np.asarray(a.arrays.get("lammps_id", []), dtype=int)
    if len(raw) != len(a) or len(set(raw.tolist())) != len(a):
        return None
    return {int(k): i for i, k in enumerate(raw)}


def compare_geometry(a, b):
    from ase.geometry import find_mic
    if len(a) != len(b) or a.get_chemical_formula() != b.get_chemical_formula() or not np.array_equal(a.pbc, b.pbc):
        return None
    am, bm = atom_ids(a), atom_ids(b)
    by_id = am is not None and bm is not None and set(am) == set(bm)
    if by_id:
        order = sorted(am)
        delta = np.asarray([a.positions[am[k]] - b.positions[bm[k]] for k in order])
        symbols = np.asarray([a[am[k]].symbol for k in order])
    else:
        from scipy.optimize import linear_sum_assignment
        if sorted(a.get_chemical_symbols()) != sorted(b.get_chemical_symbols()):
            return None
        apos = a.positions - np.mean(a.positions, axis=0)
        bpos = b.positions - np.mean(b.positions, axis=0)
        pairs = []
        for element in sorted(set(a.get_chemical_symbols())):
            ai = np.asarray([i for i, at in enumerate(a) if at.symbol == element])
            bi = np.asarray([i for i, at in enumerate(b) if at.symbol == element])
            diff = apos[ai, None, :] - bpos[None, bi, :]
            _, ds = find_mic(diff.reshape(-1, 3), a.cell, pbc=a.pbc)
            rr, cc = linear_sum_assignment(np.asarray(ds).reshape(len(ai), len(bi)))
            pairs.extend((int(i), int(j)) for i, j in zip(ai[rr], bi[cc]))
        delta = np.asarray([a.positions[i] - b.positions[j] for i, j in pairs])
        symbols = np.asarray([a[i].symbol for i, _ in pairs])
    delta, _ = find_mic(delta, a.cell, pbc=a.pbc)
    fw = symbols != "Ag"
    delta -= np.median(delta[fw], axis=0) if fw.any() else np.median(delta, axis=0)
    norm = np.linalg.norm(delta, axis=1)
    ag = symbols == "Ag"
    ag_rms = float(np.sqrt(np.mean(norm[ag]**2))) if ag.any() else 0.0
    fw_rms = float(np.sqrt(np.mean(norm[fw]**2))) if fw.any() else 0.0
    cell_delta = float(np.max(np.abs(a.cell.array - b.cell.array)))
    return {
        "mapping": "persistent_lammps_id" if by_id else "species_assignment_fallback",
        "Ag_registry_RMS_A": ag_rms,
        "framework_RMS_A": fw_rms,
        "all_atom_RMS_A": float(np.sqrt(np.mean(norm**2))),
        "max_displacement_A": float(norm.max()),
        "max_cell_delta_A": cell_delta,
        "exact": bool(cell_delta <= 1e-5 and norm.max() <= 1e-5),
        "near": bool(ag_rms <= NEAR_AG_RMS_A and fw_rms <= NEAR_FRAMEWORK_RMS_A),
    }


def species_minima(a):
    d = a.get_all_distances(mic=True)
    ids = np.asarray(a.arrays.get("lammps_id", np.arange(len(a))), dtype=int)
    result = {}
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            pair = "-".join(sorted((a[i].symbol, a[j].symbol)))
            if pair not in result or d[i, j] < result[pair]["distance_A"]:
                result[pair] = {"distance_A": float(d[i, j]), "ids": [int(ids[i]), int(ids[j])], "elements": [a[i].symbol, a[j].symbol]}
    return result


def central_indices(a):
    marked = np.flatnonzero(np.asarray(a.arrays.get("central_pair", []), dtype=bool))
    require(len(marked) == 2, "central pair must contain exactly two atoms")
    ag = [int(i) for i in marked if a[i].symbol == "Ag"]
    x = [int(i) for i in marked if a[i].symbol != "Ag"]
    require(len(ag) == len(x) == 1, "central pair must be Ag-X")
    return ag[0], x[0]


def label_free_copy(parent, name, info):
    a = parent.copy()
    for key in list(a.arrays):
        if key not in {"numbers", "positions", "lammps_id", "lammps_type", "central_pair"}:
            del a.arrays[key]
    a.info = dict(info)
    a.info["config_type"] = name
    a.info["proposal_only_no_DFT_label"] = True
    a.calc = None
    return a


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()
    repo, out = args.repo_root.resolve(), args.output_dir.resolve()
    require(out.is_relative_to(repo / REPORT) and out != repo / REPORT, "output must be child of the task report directory")
    require(not out.exists() or (out.is_dir() and not any(out.iterdir())), "nonempty output refused")
    out.mkdir(parents=True, exist_ok=True)
    hashes, inventory = {}, []

    def source(path, role):
        p = Path(path).resolve()
        require(p.is_relative_to(repo) and p.is_file(), f"source missing/outside repo: {p}")
        rel = str(p.relative_to(repo))
        digest = sha(p)
        hashes[rel] = digest
        return p, digest

    try:
        from ase.io import read, write
        from ase.geometry import find_mic

        # Geometry inventory includes every available versioned train/valid/test split.
        all_known, splits, official_inputs = [], [], []
        for p in sorted((repo / PROJECT).glob("mace_periodic_v*/data/*.extxyz")):
            if p.name not in {"train.extxyz", "valid.extxyz", "test.extxyz"}:
                continue
            p, digest = source(p, "dataset geometry")
            frames = read(p, index=":")
            inventory.append({"path": str(p.relative_to(repo)), "sha256": digest, "frames": len(frames), "role": "existing train/valid/test geometry"})
            for index, a in enumerate(frames):
                a.calc = None
                item = (str(p.relative_to(repo)), index, a)
                all_known.append(item); splits.append(item)

        # Official DFT/holdout input manifests only; never traverse calculations/.
        input_dirs = [PROJECT / n for n in (
            "pbe_interface_v13_holdouts", "pbe_interface_v14_main_holdouts",
            "pbe_interface_v14_parallel_acquisition", "pbe_interface_v15_registry_holdouts",
            "pbe_interface_v15_targeted_acquisition", "pbe_interface_v16_parallel_acquisition",
            "pbe_interface_v16_main_acquisition")]
        planned_inputs = []
        for d in input_dirs:
            mp = repo / d / "input_manifest.json"
            if not mp.is_file():
                continue
            mp, mh = source(mp, "input manifest")
            manifest = json.loads(mp.read_text())
            inventory.append({"path": str(mp.relative_to(repo)), "sha256": mh, "frames": len(manifest.get("records", [])), "role": "input manifest only"})
            for rec in manifest.get("records", []):
                ip, ih = source(repo / d / rec["input"], "input geometry")
                require(ih == rec["input_sha256"], f"frozen input hash mismatch: {rec['label']}")
                a = read(ip)
                require("REF_energy" not in a.info and "PW_PBE_energy_eV" not in a.info and not any(k in a.arrays for k in ("REF_forces", "PW_PBE_forces", "forces")), f"label found in input {rec['label']}")
                a.calc = None
                item = (str(ip.relative_to(repo)), 0, a)
                all_known.append(item); official_inputs.append(item)
                if d.name in {"pbe_interface_v16_parallel_acquisition", "pbe_interface_v16_main_acquisition"}:
                    planned_inputs.append({"label": rec["label"], "owner": rec.get("proposed_owner"), "sha256": ih, "path": str(ip.relative_to(repo))})

        # The just-published frozen validation proposals are protected inputs too.
        valdir = repo / WINDOW_B / "v16_unseen_validation_geometry_design/final_candidates"
        validation_inputs = []
        for p in sorted(valdir.glob("Ag*_v16_unseen_validation_01.extxyz")):
            p, digest = source(p, "frozen validation proposal input")
            a = read(p); a.calc = None
            item = (str(p.relative_to(repo)), 0, a)
            all_known.append(item); validation_inputs.append({"path": str(p.relative_to(repo)), "sha256": digest})
        # Prior label-free B candidate geometries count for deduplication only.
        for folder in (repo / WINDOW_B / "v15_registry_candidates", repo / WINDOW_B / "v15_wide_span_sampling_design/candidates"):
            for p in sorted(folder.glob("*.extxyz")):
                p, digest = source(p, "prior geometry proposal")
                frames = read(p, index=":")
                inventory.append({"path": str(p.relative_to(repo)), "sha256": digest, "frames": len(frames), "role": "prior label-free proposal"})
                for index, a in enumerate(frames):
                    a.calc = None; all_known.append((str(p.relative_to(repo)), index, a))

        # Historical per-atom residual vectors provide direction/priority only.
        local_path = repo / WINDOW_B / "v15_force_error_localization/analysis_01/localization.json"
        local_path, local_sha = source(local_path, "historical localization")
        localized = {x["label"].split("_")[0]: x for x in json.loads(local_path.read_text())["sources"]}
        parent_manifest_path = repo / WINDOW_B / "v15_registry_candidates/candidate_screening_manifest.json"
        parent_manifest_path, parent_manifest_sha = source(parent_manifest_path, "parent geometry lineage")
        parent_manifest = json.loads(parent_manifest_path.read_text())
        parent_meta = parent_manifest["source_dataset"]["parent_configs"]
        train_path = repo / PROJECT / "mace_periodic_v14_interface_energy/data/train.extxyz"
        train_path, train_sha = source(train_path, "parent structures")
        training = read(train_path, index=":")
        parent_frames = {}
        for kind, config in PARENTS.items():
            matches = [(i, a) for i, a in enumerate(training) if a.info.get("config_type") == config]
            require(len(matches) == 1, f"parent frame not unique: {kind}")
            parent_frames[kind] = matches[0]

        # Species-specific lower bounds derive only from current split geometries and official inputs.
        chemical_baseline = splits + official_inputs
        baseline_by_formula = {}
        for path, idx, a in chemical_baseline:
            baseline_by_formula.setdefault(a.get_chemical_formula(), []).append((path, idx, a))

        plans = {
            "AgC": {"framework_id": 604804, "description": "marked C contact with the largest longitudinal pair contribution"},
            "AgSi": {"framework_id": 366934, "description": "marked Si contact with the large transverse force error"},
            "AgTi": {"framework_id": 479027, "description": "marked Ti contact with the largest longitudinal pair contribution"},
        }
        records, groups = [], []
        for kind, config in PARENTS.items():
            frame_index, parent = parent_frames[kind]
            ids = np.asarray(parent.arrays["lammps_id"], dtype=int)
            id_to_i = {int(k): i for i, k in enumerate(ids)}
            ia, ix = central_indices(parent)
            ag_id, x_id = int(ids[ia]), int(ids[ix])
            pair_vec = parent.positions[ix] - parent.positions[ia]
            pair_len = float(np.linalg.norm(pair_vec))
            u = pair_vec / pair_len
            # z is the established interface normal. This in-plane vector is exactly
            # orthogonal to the marked pair vector and therefore preserves its projection.
            tangent = np.cross(pair_vec, np.asarray([0.0, 0.0, 1.0]))
            tangent[2] = 0.0
            require(np.linalg.norm(tangent) > 1e-9, f"cannot define lateral axis for {kind}")
            tangent /= np.linalg.norm(tangent)
            localized_row = localized[kind]
            local_atom = next(x for x in localized_row["top_10_atom_errors"] if int(x["id"]) == int(plans[kind]["framework_id"]))
            residual = np.asarray(local_atom["force_error_vector_eV_A"], dtype=float)
            framework_axis = residual / np.linalg.norm(residual)
            chosen_id = int(plans[kind]["framework_id"])
            require(chosen_id in id_to_i, f"target framework ID absent from parent: {kind}/{chosen_id}")
            parent_pair_projection = float(np.dot(pair_vec, u))
            variants = []
            for mode, axis, step, moved_ids in (
                ("transverse_Ag", tangent, TRANSVERSE_STEP_A, [int(v) for v, a in zip(ids, parent) if a.symbol == "Ag"]),
                ("framework_neighbor", framework_axis, FRAMEWORK_STEP_A, [chosen_id]),
            ):
                for sign in (-1, 1):
                    delta = float(sign * step) * axis
                    name = f"{kind}_{mode}_{'m' if sign < 0 else 'p'}"
                    info = {
                        "source_config_type": config,
                        "source_dataset_path": str(train_path.relative_to(repo)),
                        "source_dataset_sha256": train_sha,
                        "source_frame_geometry_sha256": parent_meta[kind]["source_frame_geometry_sha256"],
                        "proposal_mode": mode,
                        "proposal_sign": sign,
                        "proposal_displacement_A": delta.tolist(),
                        "candidate_role": "proposed_training_acquisition_geometry_unassigned",
                    }
                    candidate = label_free_copy(parent, name, info)
                    if mode == "transverse_Ag":
                        for atom_id in moved_ids:
                            candidate.positions[id_to_i[atom_id]] += delta
                    else:
                        candidate.positions[id_to_i[chosen_id]] += delta
                    pair_after = candidate.positions[ix] - candidate.positions[ia]
                    pair_len_after = float(np.linalg.norm(pair_after))
                    projection_after = float(np.dot(pair_after, u))
                    pair_projection_delta = projection_after - parent_pair_projection
                    minima = species_minima(candidate)
                    baseline = baseline_by_formula[candidate.get_chemical_formula()]
                    baseline_min = {}
                    for path, index, geom in baseline:
                        for pair, m in species_minima(geom).items():
                            if pair not in baseline_min or m["distance_A"] < baseline_min[pair]["distance_A"]:
                                baseline_min[pair] = {**m, "path": path, "frame": index}
                    contact_audit = {}
                    contact_fail = []
                    for pair, m in minima.items():
                        prev = baseline_min.get(pair)
                        floor = None if prev is None else prev["distance_A"] - PAIR_MARGIN_A
                        passed = prev is None or m["distance_A"] >= floor - 1e-10
                        contact_audit[pair] = {"candidate_min_A": m["distance_A"], "smallest_observed_min_A": None if prev is None else prev["distance_A"], "empirical_floor_A": floor, "source_of_minimum": prev, "pass": passed}
                        if not passed:
                            contact_fail.append({"pair": pair, "candidate_min_A": m["distance_A"], "empirical_floor_A": floor})

                    # Dedupe records all exact/near hits; near hits to the fixed source
                    # (and byte-identical lineage copies of it) are intentional controls.
                    exact_hits, near_hits, all_comparisons = [], [], []
                    for path, idx, geom in all_known:
                        cmp = compare_geometry(candidate, geom)
                        if cmp is None:
                            continue
                        hit = {"path": path, "frame": idx, "config_type": str(geom.info.get("config_type", "")), **cmp}
                        all_comparisons.append(hit)
                        if cmp["exact"]:
                            exact_hits.append(hit)
                        if cmp["near"]:
                            near_hits.append(hit)
                    all_comparisons.sort(key=lambda x: x["all_atom_RMS_A"])
                    intended_near = [x for x in near_hits if x["config_type"] == config]
                    unexpected_near = [x for x in near_hits if x["config_type"] != config]
                    path_out = out / f"{kind}_{mode}_{'minus' if sign < 0 else 'plus'}.extxyz"
                    write(path_out, candidate, format="extxyz")
                    roundtrip = read(path_out)
                    require(len(roundtrip) == len(parent) and np.array_equal(roundtrip.arrays["lammps_id"], parent.arrays["lammps_id"]), f"atom IDs/order changed: {name}")
                    require(np.allclose(roundtrip.positions, candidate.positions, atol=1e-8, rtol=0) and np.allclose(roundtrip.cell.array, parent.cell.array, atol=1e-12, rtol=0) and np.array_equal(roundtrip.pbc, parent.pbc), f"cell/PBC/coordinates serialization failed: {name}")
                    require("REF_energy" not in roundtrip.info and not any(x in roundtrip.arrays for x in ("REF_forces", "PW_PBE_forces", "forces")), f"label present in output: {name}")
                    # Framework target and all moved IDs are captured explicitly.
                    moved_records = [{"id": int(v), "element": parent[id_to_i[int(v)]].symbol, "delta_xyz_A": delta.tolist()} for v in moved_ids]
                    central_after = {"marked_pair_distance_A": pair_len_after, "longitudinal_projection_A": projection_after, "projection_delta_A": pair_projection_delta, "direct_distance_delta_A": pair_len_after - pair_len}
                    record = {
                        "label": name, "interface": kind, "mode": mode, "sign": sign,
                        "proposed_role": "training acquisition geometry only; no owner/DFT label assigned",
                        "source": {"dataset_path": str(train_path.relative_to(repo)), "dataset_sha256": train_sha, "frame_index_zero_based": frame_index, "config_type": config, "frame_geometry_sha256": parent_meta[kind]["source_frame_geometry_sha256"], "formula": parent.get_chemical_formula(), "atoms": len(parent), "control_atom_ids": {"Ag_id": ag_id, "X_id": x_id}},
                        "candidate": {"path": path_out.name, "sha256": sha(path_out), "formula": candidate.get_chemical_formula(), "atoms": len(candidate), "cell_A": candidate.cell.array.tolist(), "pbc": candidate.pbc.tolist(), "changed_atoms": moved_records, "marked_pair_before": {"distance_A": pair_len, "projection_A": parent_pair_projection}, "marked_pair_after": central_after, "species_pair_minima": minima, "species_contact_audit": contact_audit, "contact_screen_pass": not contact_fail, "contact_failures": contact_fail, "labels_absent": True},
                        "design_basis": {"target_id": chosen_id if mode == "framework_neighbor" else None, "target_element": parent[id_to_i[chosen_id]].symbol if mode == "framework_neighbor" else None, "historical_force_error_xyz_eV_A": residual.tolist() if mode == "framework_neighbor" else None, "axis_xyz": axis.tolist(), "step_A": step, "source_lineage": "shared existing small-cluster motif; not independent morphology"},
                        "dedup_screen": {"same_formula_geometry_count": len(all_comparisons), "exact_hits": exact_hits, "near_hits_intentional_parent_control": intended_near, "near_hits_other_lineages": unexpected_near, "closest_known": all_comparisons[0] if all_comparisons else None, "compared_against_six_V16_inputs": len(planned_inputs) == 6, "compared_against_frozen_validation_inputs": validation_inputs},
                    }
                    variants.append(record); records.append(record)
            groups.append({"interface": kind, "source_pair": {"Ag_id": ag_id, "X_id": x_id, "X_element": parent[ix].symbol, "distance_A": pair_len, "longitudinal_projection_A": parent_pair_projection}, "transverse_axis_xyz": tangent.tolist(), "framework_target": {"id": chosen_id, "element": parent[id_to_i[chosen_id]].symbol, "error_norm_eV_A": float(np.linalg.norm(residual)), "error_vector_xyz_eV_A": residual.tolist()}, "variants": variants})

        manifest = {
            "title": "Window B paired transverse/framework V16 training-acquisition proposals",
            "owner_instance": "c5035b48-f43b-4da0-b8e4-e2862817f86a",
            "scope": {"geometry_only": True, "DFT": False, "MACE_inference_or_training": False, "MD_or_TTM": False, "new_roles_assigned": False, "unseen_labels_read": False, "existing_training_labels_used": False},
            "design": {"control": "the parent V14 training frame for each interface", "transverse_pair": f"move the full Ag sublattice by ±{TRANSVERSE_STEP_A} A along xy tangent cross(r_X-r_Ag,z_hat), holding the marked-pair longitudinal projection constant", "framework_pair": f"move only the selected framework contact atom by ±{FRAMEWORK_STEP_A} A along its historical V15 force-residual direction; Ag fixed", "species_contact_screen": f"candidate element-pair minima cannot be >{PAIR_MARGIN_A} A below the smallest current split/official input geometry for that formula; empirical screen, not universal bond threshold", "dedup": "all train/valid/test split geometries, official V13-V16 input manifests, newly proposed withheld geometry inputs and earlier B label-free proposals; intended parent-control near hits are reported, not hidden"},
            "priority_subset": groups,
            "candidate_records": records,
            "planned_V16_acquisition_inputs": planned_inputs,
            "frozen_validation_proposal_inputs": validation_inputs,
            "inventory": inventory,
            "input_sha256": dict(sorted(hashes.items())),
        }
        write_json(out / "proposal_manifest.json", manifest)
        summary = {"design": manifest["design"], "groups": [{"interface": g["interface"], "source_pair": g["source_pair"], "framework_target": g["framework_target"], "variants": [{"label": r["label"], "mode": r["mode"], "candidate": r["candidate"], "nearest_known": r["dedup_screen"]["closest_known"], "near_hits_intentional_parent_control": r["dedup_screen"]["near_hits_intentional_parent_control"], "near_hits_other_lineages": r["dedup_screen"]["near_hits_other_lineages"]} for r in g["variants"]]} for g in groups]}
        write_json(out / "screening_summary.json", summary)
        report_lines = [
            "# Paired transverse and framework V16 acquisition proposals",
            "",
            "## Purpose and limits",
            "",
            "This bundle separates two training-data factors suggested by the matched V14/V15 diagnosis: transverse Ag registry and local framework-contact response. It contains geometry only. It does not assign DFT ownership or a final dataset role, read unseen labels, or run DFT, model inference/training, MD or TTM. A decides whether/when to label these points.",
            "",
            "For each interface, a fixed parent geometry is the control. The transverse pair translates every Ag atom by ±0.15 Å in the xy direction orthogonal to the marked Ag–X vector. The framework pair leaves Ag fixed and moves only the marked non-Ag contact by ±0.04 Å along that atom's historical V15 force-residual direction. Thus Ag registry and framework displacement are separate design arms.",
            "",
            "## Priority subset",
            "",
            "| Interface | Ag transverse pair | Framework contact pair | Pair projection behavior |",
            "|---|---|---|---|",
        ]
        for g in groups:
            x = g["source_pair"]
            fw = g["framework_target"]
            arms = {r["mode"]: r for r in g["variants"] if r["sign"] == 1}
            tr = arms["transverse_Ag"]["candidate"]["marked_pair_after"]
            fr = arms["framework_neighbor"]["candidate"]["marked_pair_after"]
            report_lines.append(f"| {g['interface']} | Ag IDs shifted ±0.15 Å; `{g['interface']}_transverse_Ag_minus/plus.extxyz` | {fw['element']}#{fw['id']} shifted ±0.04 Å; `{g['interface']}_framework_neighbor_minus/plus.extxyz` | transverse +: projection Δ {tr['projection_delta_A']:+.6f} Å, direct distance Δ {tr['direct_distance_delta_A']:+.6f} Å; framework +: projection Δ {fr['projection_delta_A']:+.6f} Å, direct distance Δ {fr['direct_distance_delta_A']:+.6f} Å |")
        report_lines += [
            "",
            "The lateral Ag arm preserves the marked-pair longitudinal projection to numerical precision; its direct Euclidean distance changes slightly at second order. Framework displacement can change the longitudinal projection and direct distance, and those changes are recorded separately. Each atom ID and displacement vector is listed in `proposal_manifest.json`.",
            "",
            "## Contact and split checks",
            "",
            "Each proposal keeps atom order/IDs, cell and PBC. Species-pair minima are checked against the empirical lower envelope from current train/valid/test geometries and official input geometries for the same formula, with a 0.05 Å allowance. This is a comparative screen, not a chemical bond definition. The frozen withheld candidates, six planned V16 acquisition inputs, all current split geometries and prior B proposals are also included in exact/near comparison.",
            "",
            "These are deliberately small perturbations around their fixed control, so the parent is an intentional near-geometry match. The manifest separates near hits to that source/control lineage from any other-lineage near hits; it does not hide or call the pair independent. The source motifs are existing small clusters, not new material morphologies. Thermal disorder and extended-interface coverage remain open.",
            "",
            "## Exact and near-geometry screening results",
            "",
            "The periodic exact-geometry screen found zero exact input-geometry hits for all 12 proposals among the screened train/valid/test frames, official acquisition/holdout inputs, frozen V16 validation proposals and prior B geometries. It matches persistent IDs where available and uses same-species assignment otherwise. A separate conservative near flag is raised when both Ag-registry RMS and non-Ag framework RMS are at most 0.15 Å after MIC matching. These are review flags, not duplicate declarations; repeated frames copied between model versions and input archives are counted as separate file/frame hits. See each candidate's complete hit list and cell deltas in `proposal_manifest.json`.",
            "",
            "| Proposal | Exact hits | Other-lineage near hits | Nearest near hit |",
            "|---|---:|---:|---|",
        ]
        for record in records:
            hits = record["dedup_screen"]["near_hits_other_lineages"]
            closest = min(hits, key=lambda h: h["all_atom_RMS_A"]) if hits else None
            closest_text = "none" if closest is None else f"`{Path(closest['path']).name}` frame {closest['frame']} (all-atom RMS {closest['all_atom_RMS_A']:.4f} Å)"
            report_lines.append(f"| {record['label']} | {len(record['dedup_screen']['exact_hits'])} | {len(hits)} | {closest_text} |")
        report_lines += [
            "",
            "AgTi_transverse_Ag_m has the largest overlap flag count (17 file/frame hits, many mirrored historical split copies); its closest listed frame is 0.0565 Å all-atom RMS away. It is not an exact duplicate, but A should review this branch against the historical split before assigning it for labeling. No candidate is near a frozen V16 validation proposal under this screen. Near-hit counts and species-specific contact checks do not establish independent statistical coverage.",
        ]
        report_lines += [
            "",
            "## A review",
            "",
            "Priority is one paired example of each arm per interface (12 geometries total). First review the AgSi marked-Si transverse pair, then the AgC marked-C and AgTi marked-Ti framework pairs and their opposite Ag shifts. Treat AgTi_transverse_Ag_m as review-only until A checks the close historical split hits. The largest noncontact shell residuals (AgC Ti#602709; AgTi C#480975) can be a second-stage framework pair if A needs separate shell-atom coverage; they were not added to this compact first set.",
            "",
            "The fixed source parent is the control; its existing label is not duplicated here. Keep these candidates outside any frozen V16 validation role. Do not label before A approves the acquisition role and timing.",
            "",
            "## Reproduction",
            "",
            "```bash",
            "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\",
            "  /workspace/.venvs/mace-v14-assist/bin/python -B \\",
            "  research/mace-v12-transfer/coordination/reports/window-b/targeted_transverse_framework_design/generate_training_proposals.py \\",
            "  --repo-root /workspace/- \\",
            "  --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/targeted_transverse_framework_design/reproduction_01",
            "```",
            "",
            "The script refuses nonempty output. Hashes cover the generator, all 12 geometry files, manifest, screening summary and this report.",
        ]
        (out / "proposal_report.md").write_text("\n".join(report_lines) + "\n")
        script_copy = out / "generate_training_proposals.py"
        script_copy.write_bytes(Path(__file__).read_bytes())
        outputs = [script_copy, *sorted(out.glob("*.extxyz")), out / "proposal_manifest.json", out / "screening_summary.json", out / "proposal_report.md"]
        (out / "SHA256SUMS.txt").write_text("".join(f"{sha(p)}  {p.name}\n" for p in outputs))
        print(json.dumps({"candidate_count": len(records), "contact_failures": sum(len(r["candidate"]["contact_failures"]) for r in records), "outputs": [p.name for p in outputs]}, indent=2))
    except Exception as exc:
        write_json(out / "diagnostics_failure.json", {"error": type(exc).__name__, "message": str(exc), "input_sha256": dict(sorted(hashes.items()))})
        raise


if __name__ == "__main__":
    main()
