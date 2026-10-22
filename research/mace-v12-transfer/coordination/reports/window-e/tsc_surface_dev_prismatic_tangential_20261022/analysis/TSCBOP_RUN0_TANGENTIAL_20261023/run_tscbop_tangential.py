#!/usr/bin/env python3
from __future__ import annotations
import csv, hashlib, importlib.util, json, os, subprocess, time
from pathlib import Path
import numpy as np
from ase.io import read

PKG = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
OLD = PKG.parent / "tsc_surface_dev_prismatic_11m20_20261021"
OLD_PAIR = OLD / "analysis/TSCBOP_RUN0_PAIR_20261022"
NEW_ID = "TSC120_DIAG_PRISM_11m20_COD9009647_10L_RESP_TIX020"
REF_ID = "TSC120_DEV_PRISM_11m20_COD9009647_10L_REF"
ATTEMPT = "attempt-0001"
REF_ATTEMPT = "attempt-0003"
LMP = Path("/workspace/.venvs/lammps-run0/bin/lmp")
RUNNER = PKG / "run_point.py"

spec = importlib.util.spec_from_file_location("old_pair", OLD_PAIR / "run_pair.py")
old_pair = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old_pair)
old_pair.OUT = OUT


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def metrics(err: np.ndarray, ids: np.ndarray, mask: np.ndarray) -> dict:
    e = err[mask]
    norms = np.linalg.norm(e, axis=1)
    selected_ids = ids[mask]
    return {
        "n_atoms": int(mask.sum()),
        "component_rmse_eV_A": float(np.sqrt(np.mean(e ** 2))),
        "vector_rmse_eV_A": float(np.sqrt(np.mean(np.sum(e ** 2, axis=1)))),
        "p95_vector_error_eV_A": float(np.percentile(norms, 95)),
        "max_vector_error_eV_A": float(norms.max()),
        "max_error_source_atom_id": int(selected_ids[np.argmax(norms)]),
    }

# Bind the new DFT result to its frozen geometry/input.
new_input = PKG / "inputs" / f"{NEW_ID}.extxyz"
new_dir = PKG / "results" / NEW_ID / ATTEMPT
new_result_path = new_dir / "RESULT.json"
new_status = json.loads((new_dir / "RUN_STATUS.json").read_text())
new_result = json.loads(new_result_path.read_text())
new_atoms = read(new_input, format="extxyz")
new_labelled = read(new_dir / "labelled.extxyz", format="extxyz")
if new_status.get("status") != "SUCCESS" or new_status.get("scf_converged") is not True:
    raise RuntimeError("New DFT result is not converged")
if sha(new_input) != new_result["input_sha256"] or old_pair.geometry_sha(new_atoms) != new_result["geometry_sha256"]:
    raise RuntimeError("New DFT input/geometry SHA mismatch")
for key in ("numbers", "positions", "cell", "pbc"):
    if not np.array_equal(getattr(new_atoms, key), getattr(new_labelled, key)):
        raise RuntimeError(f"DFT labelled geometry differs in {key}")
ids = np.asarray(new_atoms.arrays["source_atom_id"], dtype=np.int64)
if not np.array_equal(ids, new_labelled.arrays["source_atom_id"]):
    raise RuntimeError("New DFT source atom IDs mismatch")
new_dft = np.asarray(new_result["forces_eV_per_angstrom"], float)
if new_dft.shape != (len(new_atoms), 3) or not np.isfinite(new_dft).all():
    raise RuntimeError("Invalid DFT forces")

# Reuse only the already verified same-method frozen reference point and run-0 output.
ref_input = OLD / "inputs" / f"{REF_ID}.extxyz"
ref_dir = OLD / "results" / REF_ID / REF_ATTEMPT
ref_result = json.loads((ref_dir / "RESULT.json").read_text())
ref_status = json.loads((ref_dir / "RUN_STATUS.json").read_text())
ref_atoms = read(ref_input, format="extxyz")
if ref_status.get("status") != "SUCCESS" or ref_status.get("scf_converged") is not True:
    raise RuntimeError("Reference DFT is not converged")
if sha(ref_input) != ref_result["input_sha256"] or old_pair.geometry_sha(ref_atoms) != ref_result["geometry_sha256"]:
    raise RuntimeError("Reference DFT input/geometry SHA mismatch")
if not np.array_equal(ids, np.asarray(ref_atoms.arrays["source_atom_id"], dtype=np.int64)):
    raise RuntimeError("Reference and new source atom ID order differs")
if (new_atoms.get_chemical_symbols() != ref_atoms.get_chemical_symbols()
        or not np.array_equal(new_atoms.cell, ref_atoms.cell)
        or not np.array_equal(new_atoms.pbc, ref_atoms.pbc)):
    raise RuntimeError("Reference and new species/cell/PBC differ")
old_cfg = json.loads((OLD / "SCIENTIFIC_CONFIG.json").read_text())
new_cfg = json.loads((PKG / "SCIENTIFIC_CONFIG.json").read_text())
if old_cfg["method"] != new_cfg["method"]:
    raise RuntimeError("DFT reference and response methods are not identical")
ref_dft = np.asarray(ref_result["forces_eV_per_angstrom"], float)
if ref_dft.shape != (len(ref_atoms), 3) or not np.isfinite(ref_dft).all():
    raise RuntimeError("Invalid reference DFT forces")

# Run a single frozen-coordinate LAMMPS run 0 for the new geometry.
if not LMP.is_file():
    raise FileNotFoundError(LMP)
if sha(old_pair.POT) != old_pair.POT_SHA:
    raise RuntimeError("TSCBOP source hash mismatch")
run_dir = OUT / "runs" / "RESP_TIX020"
run_ids, symbols, wrapped = old_pair.make_lammps_input(new_atoms, run_dir)
env = os.environ.copy()
env.update({"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"})
started = time.monotonic()
proc = subprocess.run([str(LMP), "-in", "run.in", "-log", "lammps.log", "-screen", "screen.log", "-nocite"],
                      cwd=run_dir, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
(run_dir / "launcher.log").write_text(proc.stdout)
if proc.returncode != 0 or not (run_dir / "forces.dump").is_file():
    raise RuntimeError(f"LAMMPS run 0 failed, exit={proc.returncode}: {proc.stdout[-2000:]}")
lmp_rec = old_pair.parse_forces_dump(run_dir / "forces.dump")
ref_rec = old_pair.parse_forces_dump(OLD_PAIR / "runs/REF/forces.dump")
if set(lmp_rec) != set(ids.tolist()) or set(ref_rec) != set(ids.tolist()):
    raise RuntimeError("LAMMPS force dump IDs do not match source IDs")
idx_by_id = {int(aid): i for i, aid in enumerate(ids)}
new_lmp = np.zeros_like(new_dft)
ref_lmp = np.zeros_like(ref_dft)
max_coord_err = 0.0
for aid in ids:
    i = idx_by_id[int(aid)]
    elem, xyz, force = lmp_rec[int(aid)]
    ref_elem, _, force_ref = ref_rec[int(aid)]
    if elem != symbols[i] or ref_elem != symbols[i]:
        raise RuntimeError(f"Element mapping mismatch for source atom {aid}")
    coord_err = float(np.linalg.norm(xyz - wrapped[i]))
    max_coord_err = max(max_coord_err, coord_err)
    if coord_err > 1e-8:
        raise RuntimeError(f"LAMMPS coordinate mismatch at atom {aid}: {coord_err}")
    new_lmp[i] = force
    ref_lmp[i] = force_ref

# Absolute force comparison at the new point, plus differential response relative to REF.
abs_err = new_lmp - new_dft
dft_delta = new_dft - ref_dft
lmp_delta = new_lmp - ref_lmp
response_err = lmp_delta - dft_delta
symbols_arr = np.asarray(symbols)
all_mask = np.ones(len(ids), dtype=bool)
by_element = {e: metrics(abs_err, ids, symbols_arr == e) for e in ("Ti", "Si", "C")}
response_by_element = {e: metrics(response_err, ids, symbols_arr == e) for e in ("Ti", "Si", "C")}
moved = np.linalg.norm(new_atoms.positions - ref_atoms.positions, axis=1) > 1e-7
moved_ids = ids[moved].tolist()
if moved_ids != [4, 113]:
    raise RuntimeError(f"Unexpected displaced source IDs: {moved_ids}")
disp = new_atoms.positions[moved] - ref_atoms.positions[moved]
if not np.allclose(np.sort(disp[:, 0]), [-0.02, 0.02], atol=1e-8) or not np.allclose(disp[:, 1:], 0, atol=1e-10):
    raise RuntimeError("Tangential displacement does not match the frozen pair")
d1 = ref_atoms.get_distances(np.flatnonzero(moved)[0], np.arange(len(ids)), mic=True)
d2 = ref_atoms.get_distances(np.flatnonzero(moved)[1], np.arange(len(ids)), mic=True)
near = np.minimum(d1, d2) < 3.3

with (OUT / "TANGENTIAL_FORCE_COMPARISON.csv").open("w", newline="") as f:
    fields = ["source_atom_id", "element", "x_A", "y_A", "z_A", "dft_fx_eV_A", "dft_fy_eV_A", "dft_fz_eV_A", "tscbop_fx_eV_A", "tscbop_fy_eV_A", "tscbop_fz_eV_A", "error_fx_eV_A", "error_fy_eV_A", "error_fz_eV_A", "force_error_norm_eV_A"]
    w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
    w.writeheader()
    for i, aid in enumerate(ids):
        row = {"source_atom_id": int(aid), "element": symbols[i], "x_A": float(new_atoms.positions[i,0]), "y_A": float(new_atoms.positions[i,1]), "z_A": float(new_atoms.positions[i,2]), "force_error_norm_eV_A": float(np.linalg.norm(abs_err[i]))}
        row.update({f"dft_f{axis}_eV_A": float(new_dft[i,k]) for k,axis in enumerate("xyz")})
        row.update({f"tscbop_f{axis}_eV_A": float(new_lmp[i,k]) for k,axis in enumerate("xyz")})
        row.update({f"error_f{axis}_eV_A": float(abs_err[i,k]) for k,axis in enumerate("xyz")})
        w.writerow(row)
with (OUT / "TANGENTIAL_RESPONSE_DIFFERENCE.csv").open("w", newline="") as f:
    fields = ["source_atom_id", "element", "dft_delta_fx_eV_A", "dft_delta_fy_eV_A", "dft_delta_fz_eV_A", "tscbop_delta_fx_eV_A", "tscbop_delta_fy_eV_A", "tscbop_delta_fz_eV_A", "response_error_norm_eV_A", "within_3p3A_of_displaced_Ti"]
    w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
    w.writeheader()
    for i, aid in enumerate(ids):
        row = {"source_atom_id": int(aid), "element": symbols[i], "response_error_norm_eV_A": float(np.linalg.norm(response_err[i])), "within_3p3A_of_displaced_Ti": bool(near[i])}
        row.update({f"dft_delta_f{axis}_eV_A": float(dft_delta[i,k]) for k,axis in enumerate("xyz")})
        row.update({f"tscbop_delta_f{axis}_eV_A": float(lmp_delta[i,k]) for k,axis in enumerate("xyz")})
        w.writerow(row)

summary = {
    "status": "RUN0_AND_FORCE_MAPPING_PASS",
    "scope": "One frozen-coordinate LAMMPS run 0 for the tangential response geometry; no minimization or dynamics.",
    "role": "same COD 9009647 source family diagnostic only; no TRAIN/DEV/TEST role assignment",
    "structure_id": NEW_ID,
    "atom_count": len(ids),
    "composition": {e: int(np.sum(symbols_arr == e)) for e in ("Ti", "Si", "C")},
    "input_sha256": sha(new_input),
    "geometry_sha256": new_result["geometry_sha256"],
    "config_sha256": new_result["config_sha256"],
    "reference_input_sha256": sha(ref_input),
    "reference_geometry_sha256": ref_result["geometry_sha256"],
    "reference_config_sha256": ref_result["config_sha256"],
    "dft_reference_method_equal": True,
    "scf_converged": True,
    "scf_steps": new_result["scf_steps"],
    "potential": {"sha256": old_pair.POT_SHA, "pair_style": "tersoff", "element_order": ["Ti", "Si", "C"], "units": "metal: eV, Angstrom, eV/Angstrom"},
    "lammps": {"executable": str(LMP), "version": subprocess.run([str(LMP), "-h"], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True).stdout.splitlines()[0].strip(), "steps": 0, "runtime_seconds": time.monotonic()-started},
    "mapping": {"source_ids_match": True, "species_cell_pbc_match": True, "max_lammps_wrapped_coordinate_error_A": max_coord_err, "displaced_source_atom_ids": moved_ids, "paired_cartesian_displacements_A": {str(int(aid)): new_atoms.positions[i].tolist() for i,aid in enumerate(ids) if moved[i]}},
    "absolute_tscbop_vs_dft_at_response": {"all": metrics(abs_err, ids, all_mask), "by_element": by_element},
    "differential_force_response_error": {"definition": "(F_TSCBOP_tangent-F_TSCBOP_ref) - (F_DFT_tangent-F_DFT_ref)", "all": metrics(response_err, ids, all_mask), "by_element": response_by_element, "within_3p3A": metrics(response_err, ids, near)},
    "interpretation_limit": "A paired finite response to one symmetric tangential displacement in one COD 9009647 parent family; not an independent validation family, not a fit score, and not proof of broad TSCBOP accuracy."
}
(OUT / "SUMMARY.json").write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n")
print(json.dumps(summary, indent=2, sort_keys=True))
