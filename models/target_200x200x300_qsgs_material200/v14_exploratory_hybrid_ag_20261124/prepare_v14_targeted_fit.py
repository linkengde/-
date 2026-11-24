#!/workspace/venvs/gpaw/bin/python
"""Stage an isolated, narrowly targeted Ag correction to V13 POD.

No original V13 files are modified.  A complete 3-point Ag FCC EOS scan and
the reviewed pre-collapse liquid Ag DFT label are required before staging.
The validation liquid parent is removed as a whole to prevent parent leakage.
"""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
V13_FIT = ROOT / "fit_v13"
DIAG = ROOT / "v14_stage2_p0_ag_dft_diagnostics"
EOS = DIAG / "ag_fcc_eos_k4_equiv"
LIQUID = DIAG / "ag_liquid_precollapse" / "label.extxyz"
OUT = ROOT / "fit_v14_targeted_ag"
TRAIN_EOS = [EOS / "a_4.09A.extxyz", EOS / "a_3.90A.extxyz", EOS / "a_3.60A.extxyz"]
VALID_EOS = EOS / "a_3.75A.extxyz"


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def clean_copy(src: Path, dst: Path, candidate_id: str, role: str) -> dict:
    atoms = read(src)
    energy = float(atoms.get_potential_energy())
    forces = np.asarray(atoms.get_forces(), dtype=float)
    if not np.isfinite(energy) or forces.shape != (len(atoms), 3) or not np.isfinite(forces).all():
        raise RuntimeError(f"Invalid energy/forces in {src}")
    clean = Atoms(symbols=atoms.get_chemical_symbols(), positions=atoms.get_positions(),
                  cell=atoms.cell, pbc=atoms.pbc)
    clean.calc = SinglePointCalculator(clean, energy=energy, forces=forces)
    write(dst, clean, format="extxyz")
    return {"candidate_id": candidate_id, "role": role, "source": str(src),
            "source_sha256": sha(src), "staged_file": str(dst),
            "staged_sha256": sha(dst), "natoms": len(clean),
            "elements": sorted(set(clean.get_chemical_symbols())),
            "energy_eV": energy, "force_max_eV_A": float(np.linalg.norm(forces, axis=1).max())}


def main():
    eos_meta = [EOS / f"a_{a:.2f}A.json" for a in (4.09, 3.90, 3.60, 3.75)]
    for p in [LIQUID, VALID_EOS, *TRAIN_EOS, *eos_meta]:
        if not p.is_file():
            raise SystemExit(f"Required audited DFT label is missing: {p}")
    for p in eos_meta:
        row = json.loads(p.read_text())
        if row.get("status") != "complete" or not row.get("scf_converged") or not row.get("finite_energy_force_stress"):
            raise SystemExit(f"DFT label is not accepted for staging: {p}: {row.get('status')}")
    review = json.loads((DIAG / "ag_liquid_precollapse" / "candidate_review.json").read_text())
    if not review:
        raise SystemExit("Candidate review record is empty")

    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite existing V14 fit directory: {OUT}")
    for folder in (OUT / "training", OUT / "validation"):
        folder.mkdir(parents=True, exist_ok=False)

    staged = []
    for src in sorted((V13_FIT / "training").glob("*.xyz")):
        staged.append(clean_copy(src, OUT / "training" / src.name, src.stem, "training: retained V13"))
    if len(staged) != 27:
        raise SystemExit(f"Expected 27 V13 training labels, found {len(staged)}")

    liquid_meta = json.loads((DIAG / "ag_liquid_precollapse" / "label.json").read_text())
    staged.append(clean_copy(LIQUID, OUT / "training" / "V14_AG_LIQUID_PRECOLLAPSE_0060.xyz",
                             "V14_AG_LIQUID_PRECOLLAPSE_0060", "training: targeted DFT"))
    for src in TRAIN_EOS:
        a = src.name.split("_")[1].replace("A.extxyz", "")
        staged.append(clean_copy(src, OUT / "training" / f"V14_AG_FCC_EOS_A{a.replace('.', '')}.xyz",
                                 f"V14_AG_FCC_EOS_A{a}", "training: targeted DFT"))

    validation = []
    for src in sorted((V13_FIT / "validation").glob("*.xyz")):
        if src.stem == "V13_VAL_AG_LIQUID_01":
            continue  # entire parent trajectory reused for the targeted training label
        validation.append(clean_copy(src, OUT / "validation" / src.name, src.stem,
                                     "validation: retained V13, parent-disjoint"))
    validation.append(clean_copy(VALID_EOS, OUT / "validation" / "V14_AG_FCC_EOS_A375.xyz",
                                 "V14_AG_FCC_EOS_A3.75", "validation: held-out EOS interpolation"))
    if len(validation) != 6:
        raise SystemExit(f"Expected 6 parent-disjoint validation labels, found {len(validation)}")

    shutil.copy2(V13_FIT / "AgTiSiC_POD_v13.pod", OUT / "AgTiSiC_POD_v14.pod")
    (OUT / "fitpod.in").write_text("clear\nunits metal\natom_style atomic\nfitpod AgTiSiC_POD_v14.pod fitpod.data\n")
    (OUT / "fitpod.data").write_text(f"""file_format extxyz
file_extension xyz
path_to_training_data_set {OUT / 'training'}
path_to_test_data_set {OUT / 'validation'}
basename_for_output_files v14_targeted_ag
fitting_weight_energy 100.0
fitting_weight_force 1.0
fitting_weight_stress 0.0
error_analysis_for_training_data_set 1
error_analysis_for_test_data_set 1
energy_force_calculation_for_training_data_set 0
energy_force_calculation_for_test_data_set 0
fraction_training_data_set 1.0
fraction_test_data_set 1.0
randomize_training_data_set 0
randomize_test_data_set 0
fitting_regularization_parameter 0.1
save_pod_descriptors 0
compute_pod_descriptors 0
""")
    manifest = {
        "status": "staged_for_fit_after_label_audit",
        "training_count": len(staged), "validation_count": len(validation),
        "training_categories": {"retained_V13": 27, "precollapse_liquid_Ag": 1, "FCC_EOS": 3},
        "validation_categories": {"retained_V13_parent_disjoint": 5, "heldout_FCC_EOS": 1},
        "purged_parent": liquid_meta.get("parent_id"),
        "note": "This set targets the proven pure-Ag compression/short-range failure only. It does not repair Ti-C or prove composite interface validity; PW350 labels are consistent exploratory labels, not a claim of universal DFT convergence.",
        "staged_training": staged, "staged_validation": validation,
        "descriptor_sha256": sha(OUT / "AgTiSiC_POD_v14.pod"),
        "fit_data_sha256": sha(OUT / "fitpod.data"), "fit_input_sha256": sha(OUT / "fitpod.in"),
    }
    (OUT / "fit_data_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"training": len(staged), "validation": len(validation), "fit_dir": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
