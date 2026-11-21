#!/usr/bin/env python3
"""Refit MTP candidates and score full labels plus a masked Si-center holdout."""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "mtp_fit_candidate"
OUT = ROOT / "mtp_acquisition_refits"
BASE = Path("/workspace/potential-comparison/mace_fit_data_20261118")
MLP = Path("/workspace/potential-comparison/software/mlip-2-atomweighted-build/mlp")
INIT = BASE / "MTP10_AgCSiTi_v19v25_v24_untrained.mtp"
PYTHON = Path("/workspace/.venvs/gpaw-diagnosis/bin/python")
SCORE = BASE / "score_mlip_cfg.py"
FORCE_WEIGHTS = [0.05, 0.10, 0.20]
SCALE_BY_FORCE = [0.0]
ENERGY_GATE = 10.0
FORCE_GATE = 0.20
HOTSPOT_GATE = 0.50


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def cfg_frames(path: Path) -> list[dict]:
    text = path.read_text()
    frames = []
    for block in re.split(r"(?=^BEGIN_CFG\s*$)", text, flags=re.M):
        if not block.strip():
            continue
        size = re.search(r"^\s*Size\s*\n\s*(\d+)", block, re.M)
        atomdata = re.search(r"^\s*AtomData:.*\n((?:.*\n){%d})" % int(size.group(1)), block, re.M) if size else None
        energy = re.search(r"^\s*Energy\s*\n\s*([-+0-9.eE]+)", block, re.M)
        if not size or not atomdata or not energy:
            raise ValueError(f"malformed CFG frame in {path}")
        rows = []
        for line in atomdata.group(1).splitlines():
            fields = line.split()
            if len(fields) < 8:
                raise ValueError(f"malformed AtomData row in {path}: {line}")
            rows.append({"id": int(fields[0]), "type": int(fields[1]),
                         "force": [float(x) for x in fields[5:8]]})
        frames.append({"n": int(size.group(1)), "energy": float(energy.group(1)), "rows": rows})
    return frames


def run(args: list[str], log_path: Path) -> None:
    with log_path.open("w") as log:
        proc = subprocess.run(args, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"command returned {proc.returncode}; see {log_path}")


def fit_and_score(tag: str, training_cfg: Path, force_weight: float, scale_by_force: float,
                  center_weight: float | None, hotspot_refs: list[dict], summary: dict) -> dict:
    model = OUT / f"MTP10_v24acq_{tag}.mtp"
    run([str(MLP), "train", str(INIT), str(training_cfg),
         "--energy-weight=1", f"--force-weight={force_weight}",
         f"--scale-by-force={scale_by_force}", "--max-iter=1000",
         f"--valid-cfgs={DATA / 'dev.cfg'}", f"--trained-pot-name={model}"], OUT / f"train_{tag}.log")
    if not model.is_file():
        raise RuntimeError(f"MTP training did not write {model.name}")

    train_pred = OUT / f"train_{tag}_pred.cfg"
    dev_pred = OUT / f"dev_{tag}_pred.cfg"
    hotspot_pred = OUT / f"hotspot_centers_{tag}_pred.cfg"
    run([str(MLP), "calc-efs", str(model), str(DATA / "train.cfg"), str(train_pred)], OUT / f"calc_train_{tag}.log")
    run([str(MLP), "calc-efs", str(model), str(DATA / "dev.cfg"), str(dev_pred)], OUT / f"calc_dev_{tag}.log")
    run([str(MLP), "calc-efs", str(model), str(DATA / "hotspot_center_eval.cfg"), str(hotspot_pred)], OUT / f"calc_hotspot_{tag}.log")
    train_metrics_path = OUT / f"train_{tag}_metrics.json"
    dev_metrics_path = OUT / f"dev_{tag}_metrics.json"
    run([str(PYTHON), str(SCORE), str(DATA / "train.cfg"), str(train_pred), str(train_metrics_path)], OUT / f"score_train_{tag}.log")
    run([str(PYTHON), str(SCORE), str(DATA / "dev.cfg"), str(dev_pred), str(dev_metrics_path)], OUT / f"score_dev_{tag}.log")

    predicted_hotspots = cfg_frames(hotspot_pred)
    if len(predicted_hotspots) != len(hotspot_refs):
        raise ValueError("MTP hotspot prediction frame count mismatch")
    center_scores = []
    for ref, pred in zip(hotspot_refs, predicted_hotspots):
        idx = int(ref["center_atom_index_zero_based"])
        if pred["n"] != len(pred["rows"]):
            raise ValueError("MTP hotspot prediction atom count mismatch")
        force = pred["rows"][idx]["force"]
        target = ref["center_force_DFT_eV_A"]
        error = sum((float(force[k]) - float(target[k])) ** 2 for k in range(3)) ** 0.5
        center_scores.append({
            "index": ref["index"], "role": ref["role"], "config_type": ref["config_type"],
            "center_atom_id": ref["center_atom_id"], "prediction_eV_A": force,
            "reference_eV_A": target, "vector_error_eV_A": error,
        })

    train_metrics = json.loads(train_metrics_path.read_text())
    dev_metrics = json.loads(dev_metrics_path.read_text())
    holdout = next(item for item in center_scores if item["index"] == 1)
    result = {
        "tag": tag,
        "force_weight": force_weight,
        "scale_by_force": scale_by_force,
        "center_force_weight": center_weight,
        "fit_training_cfg": training_cfg.name,
        "model": model.name,
        "model_sha256": sha(model),
        "training": train_metrics,
        "frozen_development": dev_metrics,
        "hotspot_centers": center_scores,
        "v11_index_01_center_vector_error_eV_A": holdout["vector_error_eV_A"],
        "gates": {
            "training_energy_pass": train_metrics["energy_rmse_meV_atom"] <= ENERGY_GATE,
            "training_force_pass": train_metrics["force_component_rmse_eV_A"] <= FORCE_GATE,
            "development_energy_pass": dev_metrics["energy_rmse_meV_atom"] <= ENERGY_GATE,
            "development_force_pass": dev_metrics["force_component_rmse_eV_A"] <= FORCE_GATE,
            "heldout_center_01_pass": holdout["vector_error_eV_A"] <= HOTSPOT_GATE,
        },
    }
    result["all_gates_pass"] = all(result["gates"].values())
    (OUT / f"metrics_{tag}.json").write_text(json.dumps(result, indent=2) + "\n")
    summary["models"].append({
        "tag": tag,
        "model": model.name,
        "force_weight": force_weight,
        "scale_by_force": scale_by_force,
        "center_force_weight": center_weight,
        "train_energy_rmse_meV_atom": train_metrics["energy_rmse_meV_atom"],
        "train_force_component_rmse_eV_A": train_metrics["force_component_rmse_eV_A"],
        "dev_energy_rmse_meV_atom": dev_metrics["energy_rmse_meV_atom"],
        "dev_force_component_rmse_eV_A": dev_metrics["force_component_rmse_eV_A"],
        "center_01_error_eV_A": holdout["vector_error_eV_A"],
        "all_gates_pass": result["all_gates_pass"],
    })
    summary["status"] = "running"
    (OUT / "fit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return result


def main() -> None:
    manifest = json.loads((DATA / "candidate_manifest.json").read_text())
    for name, digest in manifest["files_sha256"].items():
        if sha(DATA / name) != digest:
            raise SystemExit(f"prepared MTP input hash mismatch: {name}")
    if not MLP.is_file() or not INIT.is_file():
        raise SystemExit("atom-weighted MLIP executable or registered MTP initialization is missing")
    OUT.mkdir(parents=True, exist_ok=True)
    summary = {
        "status": "running",
        "candidate_manifest_sha256": sha(DATA / "candidate_manifest.json"),
        "initial_model_sha256": sha(INIT),
        "weighted_mlip_source": "/workspace/potential-comparison/software/mlip-2-atomweighted-src",
        "weighted_mlip_random_initialization_seed": 20261120,
        "weighted_mlip_executable_sha256": sha(MLP),
        "weights": {"force": FORCE_WEIGHTS, "scale_by_force": SCALE_BY_FORCE,
                    "center_force": [100.0, 1000.0]},
        "thresholds": manifest["target_thresholds"],
        "models": [],
    }
    (OUT / "fit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    hotspot_refs = manifest["hotspot_center_evaluation"]

    # Baseline full-label MTP grid with the original uniform atom-force weights.
    for force_weight in FORCE_WEIGHTS:
        for scale_by_force in SCALE_BY_FORCE:
            scale_tag = str(scale_by_force).replace(".", "p")
            tag = f"F{int(round(force_weight * 100)):03d}_SF{scale_tag}"
            fit_and_score(tag, DATA / "train.cfg", force_weight, scale_by_force,
                          None, hotspot_refs, summary)

    # Center-only Si labels: crop energies and edge forces are explicitly masked.
    for center_weight in [100.0, 1000.0]:
        training_cfg = DATA / f"train_center_masked_w{int(center_weight)}.cfg"
        for force_weight in FORCE_WEIGHTS:
            tag = f"CENTER{int(center_weight)}_F{int(round(force_weight * 100)):03d}"
            fit_and_score(tag, training_cfg, force_weight, 0.0,
                          center_weight, hotspot_refs, summary)

    passing = [row for row in summary["models"] if row["all_gates_pass"]]
    summary["all_gates_pass_models"] = [row["model"] for row in passing]
    summary["selected_candidate"] = min(
        passing, key=lambda row: (row["dev_force_component_rmse_eV_A"],
                                  row["dev_energy_rmse_meV_atom"], row["center_01_error_eV_A"])
    ) if passing else None
    summary["status"] = "fit_complete_needs_review"
    summary["finished_utc"] = __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat()
    (OUT / "fit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
