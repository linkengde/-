"""Input and ownership checks for the fixed-geometry V20 slab method ladder."""
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
ENERGY = "GPAW native extrapolated energy, force_consistent=False; free_energy recorded separately"
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def inventory():
    expected = {}
    for line in (ROOT / "SHA256SUMS.txt").read_text().splitlines():
        digest, name = line.split("  ", 1)
        assert name not in expected
        path = ROOT / name
        assert path.resolve().is_relative_to(ROOT.resolve()) and path.is_file() and sha(path) == digest, f"Manifest hash mismatch: {name}"
        expected[name] = digest
    actual = {
        p.relative_to(ROOT).as_posix()
        for p in ROOT.rglob("*")
        if p.is_file() and p.name != "SHA256SUMS.txt" and not ({"calculations", "comparisons", "__pycache__"} & set(p.relative_to(ROOT).parts))
    }
    assert actual == set(expected), f"Input inventory mismatch: missing={sorted(actual-set(expected))}, extra={sorted(set(expected)-actual)}"


def validate(manifest, record, launch=False):
    if launch:
        assert manifest["launch_enabled"] is True and record["launch_enabled"] is True
    assert record["role"] == "numerical_convergence_control_not_training_or_validation"
    assert record["owner_task"] == "window-a"
    assert isinstance(record.get("owner_instance"), str) and record["owner_instance"]
    assert record["pbc"] == [True, True, False]
    assert record["poissonsolver"] == {"dipolelayer": "xy"}
    kpts = tuple(record["method"]["kpts"])
    sigma = float(record["method"]["smearing_eV"])
    assert (kpts, sigma) in {
        ((4, 4, 1), 0.1), ((6, 6, 1), 0.05), ((8, 8, 1), 0.1), ((6, 6, 1), 0.2)
    }, "Unregistered k-mesh/smearing pair"
    expected_method = dict(manifest["method"])
    expected_method.update(kpts=list(kpts), smearing_eV=sigma, pbc=record["pbc"], poissonsolver=record["poissonsolver"])
    assert record["method"] == expected_method
    assert record["method"]["xc"] == "PBE" and record["method"]["cutoff_eV"] == 500
    assert record["method"]["energy_keys"] == {"native": "energy (force_consistent=False)", "force_consistent": "free_energy"}
    assert manifest["energy_convention"] == ENERGY

    source = REPO / record["source_geometry_path"]
    inp = ROOT / record["input"]
    assert source.resolve().is_relative_to(REPO.resolve()) and sha(source) == record["source_geometry_sha256"]
    assert inp.resolve().is_relative_to(ROOT.resolve()) and sha(inp) == record["input_sha256"]
    a, parent = read(inp), read(source)
    assert a.pbc.tolist() == record["pbc"] == parent.pbc.tolist()
    assert a.get_chemical_symbols() == parent.get_chemical_symbols()
    assert np.array_equal(a.positions, parent.positions) and np.array_equal(a.cell.array, parent.cell.array)
    for key in ("lammps_id", "central_pair"):
        assert np.array_equal(a.arrays[key], parent.arrays[key])
    assert not any(k in a.info for k in ("REF_energy", "energy", "PW_PBE_energy_eV"))
    assert not any(k in a.arrays for k in ("forces", "REF_forces", "PW_PBE_forces"))
    assert np.isfinite(a.positions).all() and np.isfinite(a.cell.array).all()
    assert int(round(np.linalg.det(a.cell.array) * 1e8)) > 0
    assert len(a) == 52 and a.get_chemical_formula() == "C16Ag4Si8Ti24"
    assert record["method"]["GPAW_version"] == "26.7.0" and record["method"]["PAW_data_version"] == "1.2.1"
    return a


def load(label, launch=False):
    manifest = json.loads((ROOT / "input_manifest.json").read_text())
    record = next(row for row in manifest["records"] if row["label"] == label)
    if launch:
        inventory()
    validate(manifest, record, launch=launch)
    if launch:
        sys.path.insert(0, str(REPO / "research/mace-v12-transfer/coordination"))
        import sync_tasks as sync
        task = sync.load_task(record["owner_task"])
        sync.ensure_owner(record["owner_task"], task)
        assert record["owner_instance"] == sync.identity()
        assert label in task.get("labels", {})
        quota = Path("/sys/fs/cgroup/cpu.max").read_text().split()
        cores = len(os.sched_getaffinity(0))
        if int(os.environ.get("OMPI_COMM_WORLD_SIZE", "1")) > 1:
            assert int(os.environ["OMPI_COMM_WORLD_SIZE"]) == 4
            cpus = Path("/sys/fs/cgroup/cpuset.cpus.effective").read_text().strip()
            cores = sum((int(p.split("-")[1]) - int(p.split("-")[0]) + 1) if "-" in p else 1 for p in cpus.split(","))
        if quota[0] != "max":
            cores = min(cores, int(quota[0]) // int(quota[1]))
        assert cores >= 4, "Four physical-core execution allocation required"
        assert shutil.disk_usage(ROOT).free >= manifest["disk_budget"]["minimum_start_bytes"]
        import importlib.metadata as metadata
        import gpaw_data
        assert metadata.version("gpaw") == "26.7.0" and metadata.version("gpaw-data") == "1.2.1"
        for name, digest in record["method"]["PAW_setup_sha256"].items():
            assert sha(Path(gpaw_data.datapath()) / name) == digest
    return manifest, record
