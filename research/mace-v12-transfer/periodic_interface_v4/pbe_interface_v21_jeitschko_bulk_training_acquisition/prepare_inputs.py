#!/usr/bin/env python3
"""Build five source-traceable, geometry-only Ti3SiC2 training candidates."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
SOURCE_CIF = REPO / "research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_jeitschko1967_source_audit_20261027/Ti3SiC2_Jeitschko1967_fulltext_Fourier_variant.cif"
SOURCE_PAPER_SHA256 = "6765f99b3cc83daf3e48e4674dbeba109e5fffef8fa66cd53e266733ac4a4a63"
METHOD = {
    "xc": "PBE",
    "basis": "plane wave",
    "cutoff_eV": 500,
    "kpts": [1, 1, 1],
    "kpoint_label": "Gamma",
    "smearing_eV": 0.1,
    "pbc": [True, True, True],
    "poissonsolver": {},
    "GPAW_version": "26.7.0",
    "ASE_version": "3.29.0",
    "PAW_data_version": "1.2.1",
    "PAW_setup_sha256": {
        "Ag.PBE.gz": "bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97",
        "Ti.PBE.gz": "bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f",
        "Si.PBE.gz": "2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0",
        "C.PBE.gz": "fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983",
    },
    "parallel": {"mpi_ranks": 4, "sl_auto": True, "threads_per_rank": 1},
    "scf_thresholds": {"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
    "maxiter": 160,
    "mixer": {"type": "Pulay", "beta": 0.05, "nmaxold": 8, "weight": 100},
    "energy_keys": {
        "native": "get_potential_energy(force_consistent=False)",
        "finite_width_free_energy": "calc.results['free_energy']",
        "forces": "ASE calculator forces at FermiDirac width 0.1 eV",
    },
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_extxyz_17(path: Path, atoms) -> None:
    lattice = " ".join(format(float(value), ".17g") for value in atoms.cell.array.ravel())
    pbc = " ".join("T" if value else "F" for value in atoms.pbc)
    header = f'Lattice="{lattice}" Properties=species:S:1:pos:R:3:local_new_id:I:1 pbc="{pbc}" config_type={atoms.info["config_type"]}'
    rows = []
    for symbol, position, local_id in zip(atoms.get_chemical_symbols(), atoms.positions, atoms.arrays["local_new_id"]):
        rows.append(symbol + " " + " ".join(format(float(value), ".17g") for value in position) + f" {int(local_id)}")
    path.write_text(f"{len(atoms)}\n{header}\n" + "\n".join(rows) + "\n")


def main() -> None:
    assert SOURCE_CIF.is_file()
    source_hash = sha256(SOURCE_CIF)
    atoms = read(SOURCE_CIF)
    assert len(atoms) == 12 and atoms.get_chemical_formula() == "C4Si2Ti6"
    assert atoms.info["spacegroup"].no == 194
    atoms = atoms.repeat((2, 2, 1))
    atoms.pbc = True
    atoms.info.clear()
    for key in list(atoms.arrays):
        if key not in {"numbers", "positions"}:
            del atoms.arrays[key]
    atoms.new_array("local_new_id", np.arange(1, len(atoms) + 1, dtype=int))
    assert len(atoms) == 48 and atoms.get_chemical_formula() == "C16Si8Ti24"

    inputs = ROOT / "inputs"
    inputs.mkdir(exist_ok=True)
    cases = [{"label": "Ti3SiC2_JN1967_baseline_gamma_v21", "site_id": None, "delta_z_A": 0.0}]
    for element, local_id in [("Ti4f", 3), ("C4f", 9)]:
        for sign in (-1, 1):
            cases.append({
                "label": f"Ti3SiC2_JN1967_{element}_id{local_id:04d}_z_{'m' if sign < 0 else 'p'}002A_gamma_v21",
                "site_id": local_id,
                "delta_z_A": sign * 0.02,
            })

    records = []
    baseline_positions = atoms.positions.copy()
    for case in cases:
        frame = atoms.copy()
        if case["site_id"] is not None:
            index = int(case["site_id"]) - 1
            assert frame[index].symbol == ("Ti" if case["site_id"] == 3 else "C")
            frame.positions[index, 2] += case["delta_z_A"]
        frame.info.update({
            "config_type": case["label"],
            "proposal_role": "train_acquisition_candidate_pending_audit",
            "source_family": "Jeitschko_Nowotny_1967_single_crystal",
            "source_doi": "10.1007/BF00899949",
            "source_cif_sha256": source_hash,
            "source_site_estimate": "full-text Fourier estimates zTi=0.1357(15), zC=0.570(4); ideal full occupancy assumed",
            "parent_label": cases[0]["label"],
            "displaced_local_new_id": -1 if case["site_id"] is None else int(case["site_id"]),
            "delta_z_A": float(case["delta_z_A"]),
            "calculation_method": "PBE/PW500/Gamma/Fermi0.10/TTT; native and free energy recorded separately",
        })
        path = inputs / f"{case['label']}.extxyz"
        write_extxyz_17(path, frame)
        reread = read(path)
        assert reread.get_chemical_symbols() == frame.get_chemical_symbols()
        assert np.array_equal(reread.positions, frame.positions)
        assert np.array_equal(reread.cell.array, frame.cell.array)
        assert np.array_equal(reread.pbc, frame.pbc)
        assert np.array_equal(reread.arrays["local_new_id"], frame.arrays["local_new_id"])
        assert not any(k in reread.info for k in ("REF_energy", "energy", "PW_PBE_energy_eV"))
        assert not any(k in reread.arrays for k in ("forces", "REF_forces", "PW_PBE_forces"))
        records.append({
            **case,
            "input": str(path.relative_to(ROOT)),
            "input_sha256": sha256(path),
            "atom_count": len(frame),
            "formula": frame.get_chemical_formula(),
            "source_cif_sha256": source_hash,
            "source_paper_sha256": SOURCE_PAPER_SHA256,
            "method": METHOD,
            "role": "train_acquisition_candidate_pending_audit",
            "launch_enabled": True,
            "status": "queued",
        })
    for case in cases[1:]:
        index = case["site_id"] - 1
        assert np.allclose(atoms.positions, baseline_positions)
        changed = np.flatnonzero(np.any(np.abs(read(inputs / f"{case['label']}.extxyz").positions - baseline_positions) > 1e-12, axis=1))
        assert changed.tolist() == [index]
    manifest = {
        "task_id": "user_directed_Jeitschko1967_Ti3SiC2_backbone_DFT_acquisition_v21",
        "authorization": "User directed the assistant to calculate the missing data directly.",
        "launch_enabled": True,
        "training_enabled": False,
        "development_test_role_assignment": "pending",
        "source": {
            "doi": "10.1007/BF00899949",
            "source_cif": str(SOURCE_CIF.relative_to(REPO)),
            "source_cif_sha256": source_hash,
            "uploaded_paper_sha256": SOURCE_PAPER_SHA256,
        },
        "method": METHOD,
        "runtime_guard": {
            "one_job_at_a_time": True,
            "maximum_mpi_ranks": 4,
            "maximum_threads_per_rank": 1,
            "memory_cgroup_bytes": 34359738368,
            "minimum_starting_free_disk_bytes": 3221225472,
            "never_retry_an_OOM_job_without_revised_resource_plan": True,
        },
        "energy_convention": "Record native extrapolated total energy and finite-width free energy separately; record forces without relabeling energy convention.",
        "records": records,
    }
    (ROOT / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"source_cif_sha256": source_hash, "records": records}, indent=2))


if __name__ == "__main__":
    main()
