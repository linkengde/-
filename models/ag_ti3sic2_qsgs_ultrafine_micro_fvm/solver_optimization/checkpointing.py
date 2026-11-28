"""Atomic double-precision checkpoints and flushed step-level JSONL logs."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path

import numpy as np


def array_sha256(array: np.ndarray) -> str:
    data = np.ascontiguousarray(array)
    return hashlib.sha256(memoryview(data).cast("B")).hexdigest()


def _atomic_npy(path: Path, array: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            np.save(stream, np.asarray(array, dtype=np.float64), allow_pickle=False)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp_name, path)
    finally:
        if os.path.exists(temp_name):
            os.unlink(temp_name)


class CheckpointStore:
    """Checkpoint one solver state; static mask/source arrays are saved once."""

    def __init__(self, directory, sim, source_field, configuration):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.sim = sim
        self.source = np.ascontiguousarray(np.broadcast_to(
            np.asarray(source_field, dtype=np.float64), sim.shape
        ))
        self.mask_hash = array_sha256(sim.ag_fraction)
        self.source_hash = array_sha256(self.source)
        self.configuration = configuration
        _atomic_npy(self.directory / "ag_fraction.npy", sim.ag_fraction)
        _atomic_npy(self.directory / "source_field_W_m3.npy", self.source)

    def save(self, phase: str, next_dt_s: float, record: dict,
             iteration_profiles: list[dict] | None = None, progress: dict | None = None):
        material = self.sim.materials
        metadata = {
            "format_version": 1,
            "phase": phase,
            "time_s": float(self.sim.time_s),
            "next_dt_s": float(next_dt_s),
            "source_energy_J": float(self.sim.source_energy_J),
            "boundary_energy_out_J": float(self.sim.boundary_energy_out_J),
            "initial_stored_enthalpy_J": float(self.sim.initial_stored_enthalpy_J),
            "step_count": int(self.sim.step_count),
            "shape": list(self.sim.shape),
            "spacing_m": self.sim.spacing_m.tolist(),
            "mask_sha256": self.mask_hash,
            "source_sha256": self.source_hash,
            "material_configuration": {
                "Ag_mushy_width_K": float(material.width_K),
                "Ag_Tm_K": float(material.Tm_K),
                "Ag_latent_heat_J_kg": float(material.Lf_J_kg),
                "Ti3SiC2_limit_K": float(material.ti_Tmax_K),
            },
            "solver_configuration": {
                "class": type(self.sim).__name__,
                "nonlinear_rtol": float(self.sim.nonlinear_rtol),
                "linear_rtol": float(self.sim.linear_rtol),
                "max_newton": int(self.sim.max_newton),
                "max_pcg": int(self.sim.max_cg),
                "interface_resistance_m2_K_W": float(self.sim.interface_resistance_m2_K_W),
            },
            "run_configuration": self.configuration,
            "last_step_record": record,
            "last_newton_profiles": iteration_profiles or [],
            "progress": progress or {},
            "array_dtype": "float64",
        }
        state_path = self.directory / "latest_checkpoint.npz"
        fd, temp_name = tempfile.mkstemp(prefix="latest_checkpoint.", suffix=".tmp", dir=self.directory)
        try:
            with os.fdopen(fd, "wb") as stream:
                np.savez_compressed(
                    stream,
                    temperature_K=np.asarray(self.sim.temperature_K, dtype=np.float64),
                    enthalpy_density_J_m3=np.asarray(self.sim.enthalpy_density_J_m3, dtype=np.float64),
                    metadata_json=np.frombuffer(json.dumps(metadata, separators=(",", ":"),
                                                        allow_nan=False).encode("utf-8"), dtype=np.uint8),
                )
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temp_name, state_path)
        finally:
            if os.path.exists(temp_name):
                os.unlink(temp_name)
        return state_path

    def append_log(self, record: dict, iteration_profiles: list[dict], retry_reason: str | None = None):
        payload = {
            "accepted_step": record,
            "newton_iterations": iteration_profiles,
            "retry_reason": retry_reason,
            "time_s": float(self.sim.time_s),
            "next_dt_s": float(record["dt_s"]),
            "mask_sha256": self.mask_hash,
            "source_sha256": self.source_hash,
        }
        path = self.directory / "step_log.jsonl"
        with path.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(payload, separators=(",", ":"), allow_nan=False) + "\n")
            stream.flush()
            os.fsync(stream.fileno())

    @staticmethod
    def restore(path, sim, expected_source):
        """Restore fields and cumulative accounting after verifying all invariants."""
        path = Path(path)
        input_dir = path.parent
        mask_path = input_dir / "ag_fraction.npy"
        source_path = input_dir / "source_field_W_m3.npy"
        if not mask_path.is_file() or not source_path.is_file():
            raise FileNotFoundError("checkpoint is missing its hashed mask or source field")
        mask = np.load(mask_path, allow_pickle=False)
        source = np.load(source_path, allow_pickle=False)
        if array_sha256(mask) != array_sha256(sim.ag_fraction):
            raise ValueError("checkpoint phase mask does not match this geometry")
        source_expected = np.ascontiguousarray(np.broadcast_to(np.asarray(expected_source, dtype=np.float64), sim.shape))
        if array_sha256(source) != array_sha256(source_expected):
            raise ValueError("checkpoint source field does not match the requested run")
        with np.load(path, allow_pickle=False) as data:
            temperature = data["temperature_K"]
            enthalpy = data["enthalpy_density_J_m3"]
            metadata = json.loads(data["metadata_json"].tobytes().decode("utf-8"))
        if temperature.dtype != np.float64 or enthalpy.dtype != np.float64:
            raise ValueError("checkpoint arrays must be stored as float64")
        if list(temperature.shape) != list(sim.shape) or list(enthalpy.shape) != list(sim.shape):
            raise ValueError("checkpoint field shape mismatch")
        if metadata["mask_sha256"] != array_sha256(mask) or metadata["source_sha256"] != array_sha256(source):
            raise ValueError("checkpoint metadata hashes do not match saved inputs")
        if metadata["material_configuration"]["Ag_mushy_width_K"] != float(sim.materials.width_K):
            raise ValueError("checkpoint Ag mushy width differs from current material model")
        material = metadata["material_configuration"]
        if material["Ag_Tm_K"] != float(sim.materials.Tm_K):
            raise ValueError("checkpoint Ag melting point differs from current material model")
        if material["Ag_latent_heat_J_kg"] != float(sim.materials.Lf_J_kg):
            raise ValueError("checkpoint Ag latent heat differs from current material model")
        if material["Ti3SiC2_limit_K"] != float(sim.materials.ti_Tmax_K):
            raise ValueError("checkpoint Ti3SiC2 range differs from current material model")
        saved_config = metadata.get("run_configuration", {})
        if sim_config := saved_config.get("solver_parameters"):
            expected = {
                "nonlinear_rtol": float(sim.nonlinear_rtol), "linear_rtol": float(sim.linear_rtol),
                "max_newton": int(sim.max_newton), "max_pcg": int(sim.max_cg),
                "interface_resistance_m2_K_W": float(sim.interface_resistance_m2_K_W),
            }
            if sim_config != expected:
                raise ValueError("checkpoint solver tolerances/options differ from current run")
        sim.temperature_K[...] = temperature
        sim.enthalpy_density_J_m3[...] = enthalpy
        sim.time_s = float(metadata["time_s"])
        sim.source_energy_J = float(metadata["source_energy_J"])
        sim.boundary_energy_out_J = float(metadata["boundary_energy_out_J"])
        sim.initial_stored_enthalpy_J = float(metadata["initial_stored_enthalpy_J"])
        sim.step_count = int(metadata["step_count"])
        return metadata
