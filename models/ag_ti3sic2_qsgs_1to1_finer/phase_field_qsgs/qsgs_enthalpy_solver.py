"""Finite-volume single-temperature enthalpy solver on QSGS geometry.

This is a single-temperature continuum baseline. It is intentionally separate
from the repository's TTM-MD work: heat is deposited directly into the Ag
lattice enthalpy, so this module does not model electron/lattice nonequilibrium.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from thermophysical_model import AgEnthalpy, Ti3SiC2ThermalProperties
from qsgs_source_integration import integrate_qsgs_source_power_W
from qsgs_liquid_tracking import LiquidClusterTracker


DATA_DIR = Path(__file__).resolve().parent
ANGSTROM3_TO_M3 = 1.0e-30
EV_PER_PS_ANGSTROM2_TO_W_M2 = 1.602176634e-19 / 1.0e-12 / 1.0e-20
PHASE_FRACTION_ROUNDOFF_TOL = 1.0e-12
AG_NORMAL_BOILING_REFERENCE_K = 2435.0


@dataclass(frozen=True)
class ThermalState:
    temperature_K: np.ndarray
    ag_liquid_fraction: np.ndarray


class QSGSEnthalpyGrid:
    """Conservative finite-volume enthalpy method for the two limiting BCs.

    Projected grids retain the historical fast screening path. Native geometry
    assigns one cell to each original QSGS voxel, so face fluxes only cross
    actual shared faces. Optional native mesh refinement uniformly subdivides
    those existing voxels without regenerating the QSGS morphology.
    """

    def __init__(
        self,
        projected_npz: Path | None,
        boundary_case: str,
        peak_flux_eV_ps_A2: float = 5.0,
        gaussian_radius_A: float = 30.0,
        absorption_depth_A: float = 20.0,
        source_center_xy_A: tuple[float, float] | None = None,
        tile_xy: int = 1,
        source_geometry_npz: Path | None = None,
        native_geometry_npz: Path | None = None,
        mesh_refinement: int = 1,
    ) -> None:
        if boundary_case not in {"A", "B"}:
            raise ValueError("boundary_case must be A (adiabatic scaffold) or B (ideal contact)")
        if peak_flux_eV_ps_A2 < 0 or gaussian_radius_A <= 0 or absorption_depth_A <= 0:
            raise ValueError("Source flux must be nonnegative and lengths must be positive")
        if native_geometry_npz is None:
            if mesh_refinement != 1:
                raise ValueError("mesh_refinement applies only to native QSGS geometry")
            if projected_npz is None:
                raise ValueError("projected_npz is required unless native_geometry_npz is supplied")
            with np.load(projected_npz) as archive:
                self.ag_fraction = np.asarray(archive["ag_fraction"], dtype=float)
                self.sc_fraction = np.asarray(archive["skeleton_fraction"], dtype=float)
                self.size_A = np.asarray(archive["size_A"], dtype=float)
                self.spacing_A = np.asarray(archive["spacing_A"], dtype=float)
        else:
            if projected_npz is not None:
                raise ValueError("Choose either projected_npz or native_geometry_npz, not both")
            with np.load(native_geometry_npz, allow_pickle=False) as archive:
                skeleton = np.asarray(archive["skeleton"], dtype=bool)
                self.size_A = np.asarray(archive["size_A"], dtype=float)
            if not isinstance(mesh_refinement, (int, np.integer)) or mesh_refinement <= 0:
                raise ValueError("mesh_refinement must be a positive integer")
            if mesh_refinement > 1:
                skeleton = np.repeat(skeleton, int(mesh_refinement), axis=0)
                skeleton = np.repeat(skeleton, int(mesh_refinement), axis=1)
                skeleton = np.repeat(skeleton, int(mesh_refinement), axis=2)
            self.ag_fraction = (~skeleton).astype(float)
            self.sc_fraction = skeleton.astype(float)
            self.spacing_A = self.size_A / np.asarray(skeleton.shape, dtype=float)
            if source_geometry_npz is None:
                source_geometry_npz = native_geometry_npz
        self.mesh_refinement = int(mesh_refinement)
        # Overlap projection can leave machine-epsilon phase slivers at
        # nominally empty/full cells. They have negligible volume but can
        # create extreme per-cell temperatures because source and capacity
        # scale with the same near-zero fraction.
        self.ag_fraction = np.where(
            self.ag_fraction < PHASE_FRACTION_ROUNDOFF_TOL,
            0.0,
            np.where(
                self.ag_fraction > 1.0 - PHASE_FRACTION_ROUNDOFF_TOL,
                1.0,
                self.ag_fraction,
            ),
        )
        self.sc_fraction = 1.0 - self.ag_fraction
        self.source_geometry_npz = None if source_geometry_npz is None else Path(source_geometry_npz)
        source_ag_mask = None
        if self.source_geometry_npz is not None:
            with np.load(self.source_geometry_npz, allow_pickle=False) as source_geometry:
                source_ag_mask = ~np.asarray(source_geometry["skeleton"], dtype=bool)
                source_size_A = np.asarray(source_geometry["size_A"], dtype=float)
            if not np.allclose(source_size_A, self.size_A, rtol=0.0, atol=1.0e-8):
                raise ValueError("Binary QSGS source geometry size does not match projected domain")
        if not isinstance(tile_xy, (int, np.integer)) or tile_xy <= 0:
            raise ValueError("tile_xy must be a positive integer")
        self.tile_xy = int(tile_xy)
        if self.tile_xy > 1:
            self.ag_fraction = np.tile(self.ag_fraction, (self.tile_xy, self.tile_xy, 1))
            self.sc_fraction = np.tile(self.sc_fraction, (self.tile_xy, self.tile_xy, 1))
            self.size_A = self.size_A * np.asarray([self.tile_xy, self.tile_xy, 1.0])
            if source_ag_mask is not None:
                source_ag_mask = np.tile(source_ag_mask, (self.tile_xy, self.tile_xy, 1))
        if self.ag_fraction.ndim != 3 or self.sc_fraction.shape != self.ag_fraction.shape:
            raise ValueError("Projected phase fractions must be matching three-dimensional arrays")
        if not np.allclose(self.ag_fraction + self.sc_fraction, 1.0, atol=2.0e-12):
            raise ValueError("Projected Ag and scaffold fractions must sum to one")
        if np.any(self.ag_fraction < 0.0) or np.any(self.sc_fraction < 0.0):
            raise ValueError("Phase fractions cannot be negative")

        self.boundary_case = boundary_case
        self.ag = AgEnthalpy()
        self.sc = Ti3SiC2ThermalProperties()
        self.peak_flux_W_m2 = peak_flux_eV_ps_A2 * EV_PER_PS_ANGSTROM2_TO_W_M2
        self.gaussian_radius_A = float(gaussian_radius_A)
        self.absorption_depth_A = float(absorption_depth_A)
        self.cell_volume_m3 = float(np.prod(self.spacing_A) * ANGSTROM3_TO_M3)
        self.temperature_shape = self.ag_fraction.shape
        self.source_power_in_roundoff_only_cells_W = 0.0
        if source_ag_mask is None:
            x = (np.arange(self.temperature_shape[0]) + 0.5) * self.spacing_A[0]
            y = (np.arange(self.temperature_shape[1]) + 0.5) * self.spacing_A[1]
            z = (np.arange(self.temperature_shape[2]) + 0.5) * self.spacing_A[2]
            cx, cy = source_center_xy_A or (self.size_A[0] / 2.0, self.size_A[1] / 2.0)
            self.source_center_xy_A = (float(cx), float(cy))
            rr2 = (x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2
            depth = self.size_A[2] - z[None, None, :]
            spatial = np.exp(-rr2 / self.gaussian_radius_A**2) * np.exp(-depth / self.absorption_depth_A)
            self.source_W_m3 = (
                self.ag_fraction
                * self.peak_flux_W_m2
                / (self.absorption_depth_A * 1.0e-10)
                * spatial
            )
            self.actual_source_power_W = float(self.source_W_m3.sum() * self.cell_volume_m3)
            self.source_cell_power_W = self.source_W_m3 * self.cell_volume_m3
            self.source_integration_method = "cell-center source times projected Ag fraction"
        else:
            self.source_center_xy_A = (
                tuple(float(value) for value in source_center_xy_A)
                if source_center_xy_A is not None
                else (float(self.size_A[0] / 2.0), float(self.size_A[1] / 2.0))
            )
            integrated_cell_power_W = integrate_qsgs_source_power_W(
                source_ag_mask,
                self.temperature_shape,
                self.size_A,
                self.peak_flux_W_m2,
                self.gaussian_radius_A,
                self.absorption_depth_A,
                source_center_xy_A,
            )
            active_ag_cell = self.ag_fraction >= PHASE_FRACTION_ROUNDOFF_TOL
            self.source_power_in_roundoff_only_cells_W = float(
                integrated_cell_power_W[~active_ag_cell].sum()
            )
            self.source_cell_power_W = np.where(active_ag_cell, integrated_cell_power_W, 0.0)
            self.source_W_m3 = self.source_cell_power_W / self.cell_volume_m3
            self.actual_source_power_W = float(self.source_cell_power_W.sum())
            self.source_integration_method = (
                "exact separable integral over binary QSGS Ag voxels; contributions to "
                "round-off-only projected cells are omitted"
            )
        self._solid_at_tm_J_m3 = self._solid_volumetric_enthalpy(
            np.full(self.temperature_shape, self.ag.melting_temperature_K)
        )
        self._latent_volumetric_J_m3 = (
            self.ag_fraction * self.ag.reference_density_kg_m3 * self.ag.latent_heat_J_kg
        )
        self._solid_capacity_at_300_J_m3_K = (
            self.ag_fraction * self.ag.reference_density_kg_m3 * self.ag.solid_cp_J_kg_K[0]
        )
        self._solid_capacity_slope_J_m3_K2 = self.ag_fraction * self.ag.reference_density_kg_m3 * (
            (self.ag.solid_cp_J_kg_K[1] - self.ag.solid_cp_J_kg_K[0])
            / (self.ag.solid_temperature_K[1] - self.ag.solid_temperature_K[0])
        )
        if self.boundary_case == "B":
            sc_cp0 = self.sc.cp_mass_J_kg_K[0]
            sc_slope = (self.sc.cp_mass_J_kg_K[1] - sc_cp0) / (
                self.sc.temperature_K[1] - self.sc.temperature_K[0]
            )
            self._solid_capacity_at_300_J_m3_K += self.sc_fraction * self.sc.density_kg_m3 * sc_cp0
            self._solid_capacity_slope_J_m3_K2 += self.sc_fraction * self.sc.density_kg_m3 * sc_slope
        self.state = self.state_from_enthalpy_density(np.zeros(self.temperature_shape, dtype=float))
        self.enthalpy_density_J_m3 = np.zeros(self.temperature_shape, dtype=float)
        self.time_s = 0.0
        self.injected_energy_J = 0.0
        self.max_temperature_seen_K = 300.0
        self.last_run_step_limit_s: float | None = None

    @staticmethod
    def _sha256_file(path: Path) -> str:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    @staticmethod
    def _sha256_array(array: np.ndarray) -> str:
        return hashlib.sha256(np.ascontiguousarray(array).view(np.uint8)).hexdigest()

    def checkpoint_identity(self) -> dict[str, object]:
        geometry_hash = (
            self._sha256_file(self.source_geometry_npz)
            if self.source_geometry_npz is not None
            else None
        )
        table_path = DATA_DIR / "thermophysical_candidates.csv"
        config_path = DATA_DIR / "thermophysical_candidates.json"
        model_path = DATA_DIR / "thermophysical_model.py"
        return {
            "schema_version": 1,
            "boundary_case": self.boundary_case,
            "boundary_conditions": {"x": "periodic", "y": "periodic", "z": "insulated"},
            "grid_shape_xyz": list(self.temperature_shape),
            "spacing_A": self.spacing_A.tolist(),
            "domain_size_A": self.size_A.tolist(),
            "ag_fraction_sha256": self._sha256_array(self.ag_fraction),
            "geometry_sha256": geometry_hash,
            "source_geometry_sha256": geometry_hash,
            "source_integration_method": self.source_integration_method,
            "source_peak_flux_W_m2": self.peak_flux_W_m2,
            "source_gaussian_radius_A": self.gaussian_radius_A,
            "source_absorption_depth_A": self.absorption_depth_A,
            "source_center_xy_A": list(self.source_center_xy_A),
            "source_actual_power_W": self.actual_source_power_W,
            "tile_xy": self.tile_xy,
            "mesh_refinement": self.mesh_refinement,
            "thermophysical_candidates_csv_sha256": self._sha256_file(table_path),
            "thermophysical_candidates_json_sha256": self._sha256_file(config_path),
            "thermophysical_model_code_sha256": self._sha256_file(model_path),
            "enthalpy_reference_temperature_K": 300.0,
            "phase_model": "Ag solid sensible enthalpy + isothermal fusion latent heat + liquid sensible enthalpy",
        }

    def save_checkpoint(
        self,
        path: Path,
        tracker: LiquidClusterTracker | None = None,
    ) -> Path:
        """Atomically save the full restart state and its validation metadata."""
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        tracker_metadata = None if tracker is None else tracker.to_metadata()
        tracker_labels = (
            np.zeros((0, 0, 0), dtype=np.int32)
            if tracker is None or tracker.previous_track_labels is None
            else tracker.previous_track_labels.astype(np.int32, copy=False)
        )
        arrays = {
            "enthalpy_density_J_m3": np.asarray(self.enthalpy_density_J_m3, dtype=np.float64),
            "temperature_K": np.asarray(self.state.temperature_K, dtype=np.float64),
            "ag_liquid_fraction": np.asarray(self.state.ag_liquid_fraction, dtype=np.float64),
            "tracker_track_labels": tracker_labels,
        }
        total_energy = float(self.enthalpy_density_J_m3.sum() * self.cell_volume_m3)
        latent_energy = self._ag_latent_energy_J()
        energy_error = total_energy - self.injected_energy_J
        metadata = {
            "checkpoint_schema_version": 1,
            "identity": self.checkpoint_identity(),
            "identity_sha256": hashlib.sha256(
                json.dumps(self.checkpoint_identity(), sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest(),
            "time_s": self.time_s,
            "time_ps": self.time_s * 1.0e12,
            "max_temperature_seen_K": self.max_temperature_seen_K,
            "last_run_step_limit_s": self.last_run_step_limit_s,
            "injected_energy_J": self.injected_energy_J,
            "actual_source_power_W": self.actual_source_power_W,
            "total_enthalpy_J": total_energy,
            "ag_latent_energy_J": latent_energy,
            "sensible_enthalpy_J": total_energy - latent_energy,
            "energy_balance_error_J": energy_error,
            "relative_energy_balance_error": (
                energy_error / self.injected_energy_J if self.injected_energy_J != 0.0 else 0.0
            ),
            "tracker_metadata": tracker_metadata,
            "tracker_labels_present": bool(tracker is not None and tracker.previous_track_labels is not None),
            "array_sha256": {name: self._sha256_array(value) for name, value in arrays.items()},
        }
        temporary = target.with_name(target.name + ".tmp")
        try:
            with temporary.open("wb") as handle:
                np.savez_compressed(handle, **arrays, metadata_json=np.asarray(json.dumps(metadata)))
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, target)
        finally:
            if temporary.exists():
                temporary.unlink()
        return target

    def load_checkpoint(
        self,
        path: Path,
        tracker: LiquidClusterTracker | None = None,
    ) -> dict[str, object]:
        """Restore a checkpoint only when its geometry and physics match exactly."""
        with np.load(path, allow_pickle=False) as archive:
            metadata = json.loads(str(archive["metadata_json"].item()))
            if int(metadata.get("checkpoint_schema_version", -1)) != 1:
                raise ValueError("Unsupported enthalpy checkpoint schema")
            expected_identity = self.checkpoint_identity()
            expected_digest = hashlib.sha256(
                json.dumps(expected_identity, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest()
            if metadata.get("identity_sha256") != expected_digest or metadata.get("identity") != expected_identity:
                raise ValueError("Checkpoint geometry, source, boundary, or thermal-property identity does not match")
            arrays = {
                name: np.asarray(archive[name]).copy()
                for name in (
                    "enthalpy_density_J_m3",
                    "temperature_K",
                    "ag_liquid_fraction",
                    "tracker_track_labels",
                )
            }
        for name, array in arrays.items():
            if self._sha256_array(array) != metadata["array_sha256"].get(name):
                raise ValueError(f"Checkpoint array integrity check failed: {name}")
        enthalpy = arrays["enthalpy_density_J_m3"]
        temperature = arrays["temperature_K"]
        liquid = arrays["ag_liquid_fraction"]
        if enthalpy.shape != self.temperature_shape or temperature.shape != self.temperature_shape or liquid.shape != self.temperature_shape:
            raise ValueError("Checkpoint fields do not match the active grid")
        if not all(np.all(np.isfinite(array)) for array in (enthalpy, temperature, liquid)):
            raise ValueError("Checkpoint contains non-finite thermal fields")
        restored_state = self.state_from_enthalpy_density(enthalpy)
        if not np.allclose(temperature, restored_state.temperature_K, rtol=0.0, atol=1.0e-8):
            raise ValueError("Checkpoint temperature field is inconsistent with its enthalpy field")
        if not np.allclose(liquid, restored_state.ag_liquid_fraction, rtol=0.0, atol=1.0e-12):
            raise ValueError("Checkpoint liquid fraction is inconsistent with its enthalpy field")
        self.enthalpy_density_J_m3 = enthalpy
        self.state = restored_state
        self.time_s = float(metadata["time_s"])
        self.injected_energy_J = float(metadata["injected_energy_J"])
        self.max_temperature_seen_K = float(metadata["max_temperature_seen_K"])
        self.last_run_step_limit_s = (
            None
            if metadata.get("last_run_step_limit_s") is None
            else float(metadata["last_run_step_limit_s"])
        )
        total_energy = float(enthalpy.sum() * self.cell_volume_m3)
        error = total_energy - self.injected_energy_J
        stored_error = float(metadata["energy_balance_error_J"])
        if not np.isclose(error, stored_error, rtol=1.0e-10, atol=1.0e-30):
            raise ValueError("Checkpoint energy ledger is inconsistent with its enthalpy field")
        if "ag_latent_energy_J" in metadata:
            latent_energy = self._ag_latent_energy_J()
            if not np.isclose(latent_energy, float(metadata["ag_latent_energy_J"]), rtol=1.0e-12, atol=1.0e-30):
                raise ValueError("Checkpoint latent-energy ledger is inconsistent with its liquid fraction field")
        if tracker is None:
            if metadata.get("tracker_metadata") is not None:
                raise ValueError("Checkpoint contains liquid tracking state; construct a tracker to restore it")
        else:
            tracker_metadata = metadata.get("tracker_metadata")
            labels = arrays["tracker_track_labels"] if metadata.get("tracker_labels_present") else None
            if tracker_metadata is None:
                raise ValueError("Checkpoint contains no liquid tracking state")
            if labels is not None and labels.shape != self.temperature_shape:
                raise ValueError("Checkpoint tracker labels do not match the active grid")
            tracker.load_metadata(tracker_metadata, labels)
        return {
            "time_ps": self.time_s * 1.0e12,
            "injected_energy_J": self.injected_energy_J,
            "total_enthalpy_J": total_energy,
            "energy_balance_error_J": error,
        }

    @property
    def melting_temperature_K(self) -> float:
        return self.ag.melting_temperature_K

    def _solid_volumetric_enthalpy(self, temperature_K: np.ndarray) -> np.ndarray:
        """Mixture sensible enthalpy relative to 300 K, below Ag melting."""
        temp = np.asarray(temperature_K, dtype=float)
        ag_t = np.minimum(temp, self.ag.melting_temperature_K)
        h_ag = np.asarray(self.ag._integrate_solid_cp(ag_t))
        ag_part = self.ag_fraction * self.ag.reference_density_kg_m3 * h_ag
        if self.boundary_case == "A":
            return ag_part

        cp0 = self.sc.cp_mass_J_kg_K[0]
        slope = (self.sc.cp_mass_J_kg_K[1] - cp0) / (self.sc.temperature_K[1] - self.sc.temperature_K[0])
        delta = temp - self.sc.temperature_K[0]
        h_sc = cp0 * delta + 0.5 * slope * delta**2
        h_sc = np.where(temp < self.sc.temperature_K[0], cp0 * delta, h_sc)
        return ag_part + self.sc_fraction * self.sc.density_kg_m3 * h_sc

    def enthalpy_density_from_state(self, temperature_K: np.ndarray, liquid_fraction: np.ndarray) -> np.ndarray:
        temp, liquid = np.broadcast_arrays(np.asarray(temperature_K, dtype=float), np.asarray(liquid_fraction, dtype=float))
        if np.any((liquid < 0.0) | (liquid > 1.0)):
            raise ValueError("Ag liquid fraction must lie in [0,1]")
        tm = self.melting_temperature_K
        solid_at_tm = self._solid_at_tm_J_m3
        latent_vol = self._latent_volumetric_J_m3
        h_solid = self._solid_volumetric_enthalpy(temp)
        h_plateau = solid_at_tm + latent_vol * liquid
        h_liquid = solid_at_tm + latent_vol + self.ag_fraction * self.ag.reference_density_kg_m3 * float(self.ag.liquid_cp_J_kg_K[0]) * (temp - tm)
        if self.boundary_case == "B":
            dtemp = temp - tm
            slope_sc = (self.sc.cp_mass_J_kg_K[1] - self.sc.cp_mass_J_kg_K[0]) / (self.sc.temperature_K[1] - self.sc.temperature_K[0])
            cp_sc_tm = float(self.sc.cp_J_kg_K(tm))
            h_liquid += self.sc_fraction * self.sc.density_kg_m3 * (cp_sc_tm * dtemp + 0.5 * slope_sc * dtemp**2)
        out = np.where(temp < tm, h_solid, np.where(temp > tm, h_liquid, h_plateau))
        return out

    def state_from_enthalpy_density(self, enthalpy_density_J_m3: np.ndarray) -> ThermalState:
        h = np.asarray(enthalpy_density_J_m3, dtype=float)
        if h.shape != self.temperature_shape:
            h = np.broadcast_to(h, self.temperature_shape)
        temp = np.full(self.temperature_shape, 300.0, dtype=float)
        liquid = np.zeros(self.temperature_shape, dtype=float)
        ag_rho = self.ag.reference_density_kg_m3
        c0 = self._solid_capacity_at_300_J_m3_K
        slope = self._solid_capacity_slope_J_m3_K2

        below_300 = h < 0.0
        has_capacity = c0 > 0.0
        temp[below_300 & has_capacity] = 300.0 + h[below_300 & has_capacity] / c0[below_300 & has_capacity]
        solid_at_tm = self._solid_at_tm_J_m3
        solid = (h >= 0.0) & (h < solid_at_tm) & has_capacity
        target = h[solid]
        c = c0[solid]
        sl = slope[solid]
        if target.size:
            delta = target / c
            curved = np.abs(sl) > 1.0e-18
            delta[curved] = 2.0 * target[curved] / (c[curved] + np.sqrt(c[curved] ** 2 + 2.0 * sl[curved] * target[curved]))
            temp[solid] = 300.0 + delta

        latent_vol = self._latent_volumetric_J_m3
        liquidus_enthalpy = solid_at_tm + latent_vol
        plateau = (h >= solid_at_tm) & (h <= liquidus_enthalpy) & (latent_vol > 0.0)
        temp[plateau] = self.melting_temperature_K
        liquid[plateau] = np.clip((h[plateau] - solid_at_tm[plateau]) / latent_vol[plateau], 0.0, 1.0)
        no_ag_at_tm = (latent_vol <= 0.0) & (solid_at_tm > 0.0) & np.isclose(
            h, solid_at_tm, rtol=2.0e-14, atol=1.0e-8
        )
        temp[no_ag_at_tm] = self.melting_temperature_K

        above = (h > liquidus_enthalpy) & has_capacity
        if np.any(above):
            dH = h[above] - liquidus_enthalpy[above]
            cp_liquid = float(self.ag.liquid_cp_J_kg_K[0])
            c_tm = self.ag_fraction[above] * ag_rho * cp_liquid
            slope_liquid = np.zeros_like(c_tm)
            if self.boundary_case == "B":
                sc_cp_tm = float(self.sc.cp_J_kg_K(self.melting_temperature_K))
                sc_slope = (self.sc.cp_mass_J_kg_K[1] - self.sc.cp_mass_J_kg_K[0]) / (self.sc.temperature_K[1] - self.sc.temperature_K[0])
                c_tm += self.sc_fraction[above] * self.sc.density_kg_m3 * sc_cp_tm
                slope_liquid = self.sc_fraction[above] * self.sc.density_kg_m3 * sc_slope
            delta = dH / c_tm
            curved = np.abs(slope_liquid) > 1.0e-18
            delta[curved] = 2.0 * dH[curved] / (c_tm[curved] + np.sqrt(c_tm[curved] ** 2 + 2.0 * slope_liquid[curved] * dH[curved]))
            temp[above] = self.melting_temperature_K + delta
            liquid[above] = 1.0

        if self.boundary_case == "B" and np.any(temp[self.sc_fraction > 0.0] > self.sc.temperature_K[-1] + 1.0e-9):
            raise ValueError("Ti3SiC2 candidate table exceeds 1473.15 K; add high-temperature data before continuing")
        return ThermalState(temp, liquid)

    def conductivity_W_m_K(self, state: ThermalState | None = None) -> np.ndarray:
        state = self.state if state is None else state
        ag_k = np.asarray(self.ag.conductivity_W_m_K(state.temperature_K, state.ag_liquid_fraction))
        result = self.ag_fraction * ag_k
        if self.boundary_case == "B":
            sc_k = np.asarray(self.sc.conductivity_W_m_K(state.temperature_K))
            result = result + self.sc_fraction * sc_k
        return result

    def sensible_capacity_J_m3_K(self, state: ThermalState | None = None) -> np.ndarray:
        state = self.state if state is None else state
        ag_solid = np.asarray(self.ag.cp_solid_J_kg_K(np.minimum(state.temperature_K, self.melting_temperature_K)))
        ag_liquid = np.asarray(self.ag.cp_liquid_J_kg_K(np.maximum(state.temperature_K, self.melting_temperature_K)))
        ag_cp = (1.0 - state.ag_liquid_fraction) * ag_solid + state.ag_liquid_fraction * ag_liquid
        capacity = self.ag_fraction * self.ag.reference_density_kg_m3 * ag_cp
        if self.boundary_case == "B":
            capacity += self.sc_fraction * self.sc.density_kg_m3 * np.asarray(self.sc.cp_J_kg_K(state.temperature_K))
        return capacity

    def stable_timestep_s(self, safety: float = 0.45) -> float:
        """Conservative explicit diffusion limit using the current state."""
        k = self.conductivity_W_m_K()
        capacity = self.sensible_capacity_J_m3_K()
        inverse_spacing2 = 1.0 / (self.spacing_A * 1.0e-10) ** 2
        face_sum = np.zeros(self.temperature_shape, dtype=float)
        for axis in range(3):
            right_k = np.roll(k, -1, axis=axis) if axis < 2 else np.concatenate((k[:, :, 1:], np.zeros_like(k[:, :, :1])), axis=2)
            left_k = np.roll(k, 1, axis=axis) if axis < 2 else np.concatenate((np.zeros_like(k[:, :, :1]), k[:, :, :-1]), axis=2)
            right = np.divide(2.0 * k * right_k, k + right_k, out=np.zeros_like(k), where=(k + right_k) > 0.0)
            left = np.divide(2.0 * k * left_k, k + left_k, out=np.zeros_like(k), where=(k + left_k) > 0.0)
            face_sum += (right + left) * inverse_spacing2[axis]
        rate = np.divide(face_sum, capacity, out=np.zeros_like(face_sum), where=capacity > 0.0)
        max_rate = float(rate.max())
        if max_rate <= 0.0:
            return np.inf
        return safety / max_rate

    def diffusion_power_density_W_m3(self, temperature_K: np.ndarray, conductivity_W_m_K: np.ndarray) -> np.ndarray:
        """Pairwise-conservative periodic x/y, insulated z diffusion operator."""
        result = np.zeros(self.temperature_shape, dtype=float)
        for axis in (0, 1):
            dx = self.spacing_A[axis] * 1.0e-10
            t_right = np.roll(temperature_K, -1, axis=axis)
            k_right = np.roll(conductivity_W_m_K, -1, axis=axis)
            kface = np.divide(
                2.0 * conductivity_W_m_K * k_right,
                conductivity_W_m_K + k_right,
                out=np.zeros_like(conductivity_W_m_K),
                where=(conductivity_W_m_K + k_right) > 0.0,
            )
            result += kface * (t_right - temperature_K) / dx**2
            result += np.roll(kface, 1, axis=axis) * (np.roll(temperature_K, 1, axis=axis) - temperature_K) / dx**2
        dz = self.spacing_A[2] * 1.0e-10
        kface = np.divide(
            2.0 * conductivity_W_m_K[:, :, :-1] * conductivity_W_m_K[:, :, 1:],
            conductivity_W_m_K[:, :, :-1] + conductivity_W_m_K[:, :, 1:],
            out=np.zeros_like(conductivity_W_m_K[:, :, :-1]),
            where=(conductivity_W_m_K[:, :, :-1] + conductivity_W_m_K[:, :, 1:]) > 0.0,
        )
        flux_term = kface * (temperature_K[:, :, 1:] - temperature_K[:, :, :-1]) / dz**2
        result[:, :, :-1] += flux_term
        result[:, :, 1:] -= flux_term
        return result

    def advance(
        self,
        dt_s: float,
        source_on: bool = True,
        enforce_stability: bool = True,
    ) -> dict[str, float]:
        if dt_s <= 0.0 or not np.isfinite(dt_s):
            raise ValueError("dt_s must be finite and positive")
        if enforce_stability and dt_s > self.stable_timestep_s() * (1.0 + 1.0e-10):
            raise ValueError("Requested timestep exceeds the current explicit diffusion stability limit")
        k = self.conductivity_W_m_K()
        diffusion = self.diffusion_power_density_W_m3(self.state.temperature_K, k)
        source = self.source_W_m3 if source_on else 0.0
        self.enthalpy_density_J_m3 = self.enthalpy_density_J_m3 + dt_s * (diffusion + source)
        self.time_s += dt_s
        if source_on:
            self.injected_energy_J += self.actual_source_power_W * dt_s
        self.state = self.state_from_enthalpy_density(self.enthalpy_density_J_m3)
        self.max_temperature_seen_K = max(self.max_temperature_seen_K, float(self.state.temperature_K.max()))
        total_energy = float(self.enthalpy_density_J_m3.sum() * self.cell_volume_m3)
        error = total_energy - self.injected_energy_J
        return {
            "time_ps": self.time_s * 1.0e12,
            "dt_fs": dt_s * 1.0e15,
            "max_temperature_K": float(self.state.temperature_K.max()),
            "mean_ag_temperature_K": self._ag_weighted_mean_temperature(),
            "max_ag_liquid_fraction": float(self.state.ag_liquid_fraction.max()),
            "ag_liquid_volume_fraction": self._liquid_ag_volume_fraction(),
            "total_energy_J": total_energy,
            "injected_energy_J": self.injected_energy_J,
            "energy_balance_error_J": error,
        }

    def _ag_weighted_mean_temperature(self) -> float:
        denominator = float(self.ag_fraction.sum())
        if denominator <= 0.0:
            return 300.0
        return float(np.sum(self.ag_fraction * self.state.temperature_K) / denominator)

    def _liquid_ag_volume_fraction(self) -> float:
        total = float(self.ag_fraction.sum())
        return float(np.sum(self.ag_fraction * self.state.ag_liquid_fraction) / total) if total else 0.0

    def _ag_latent_energy_J(self) -> float:
        return float(
            np.sum(self._latent_volumetric_J_m3 * self.state.ag_liquid_fraction)
            * self.cell_volume_m3
        )

    def run(
        self,
        duration_ps: float,
        source_on: bool,
        sample_interval_ps: float = 0.1,
        maximum_dt_s: float | None = None,
        track_liquid_connectivity: bool = False,
        liquid_threshold: float = 0.5,
        progress: bool = False,
        tracker: LiquidClusterTracker | None = None,
        checkpoint_times_ps: tuple[float, ...] | list[float] = (),
        checkpoint_dir: Path | None = None,
        checkpoint_prefix: str = "checkpoint",
    ) -> list[dict[str, float]]:
        if duration_ps < 0.0 or sample_interval_ps <= 0.0:
            raise ValueError("Run duration must be nonnegative and sample interval positive")
        if duration_ps == 0.0 and not checkpoint_times_ps:
            return []
        if maximum_dt_s is not None and (not np.isfinite(maximum_dt_s) or maximum_dt_s <= 0.0):
            raise ValueError("maximum_dt_s must be finite and positive")
        if (track_liquid_connectivity or tracker is not None) and self.source_geometry_npz is None:
            raise ValueError("Liquid connectivity tracking requires the native QSGS geometry")
        if not 0.0 < liquid_threshold <= 1.0:
            raise ValueError("liquid_threshold must be in (0, 1]")
        if tracker is not None and not np.isclose(tracker.liquid_threshold, liquid_threshold, rtol=0.0, atol=0.0):
            raise ValueError("Tracker threshold differs from the requested liquid threshold")
        if track_liquid_connectivity and tracker is None:
            tracker = LiquidClusterTracker(
                self.spacing_A,
                self.size_A,
                liquid_threshold=liquid_threshold,
            )
        checkpoint_times = sorted(set(float(value) for value in checkpoint_times_ps))
        if any(not np.isfinite(value) or value < 0.0 for value in checkpoint_times):
            raise ValueError("Checkpoint times must be finite nonnegative global simulation times")
        if checkpoint_times and checkpoint_dir is None:
            raise ValueError("checkpoint_dir is required when checkpoint times are provided")
        end_s = self.time_s + duration_ps * 1.0e-12
        sample_interval_s = sample_interval_ps * 1.0e-12
        sample_index = int(np.floor((self.time_s + 1.0e-24) / sample_interval_s)) + 1
        next_sample = sample_index * sample_interval_s
        pending_checkpoints = [value * 1.0e-12 for value in checkpoint_times if value * 1.0e-12 > self.time_s + 1.0e-24]
        # Over the candidate range k(T)/Cp(T) decreases from its 300 K value
        # in both phases; a fixed initial stable step is therefore conservative
        # and avoids recomputing the full-grid diffusion bound every substep.
        stable_dt = self.stable_timestep_s()
        step_limit = stable_dt if maximum_dt_s is None else min(stable_dt, maximum_dt_s)
        if self.last_run_step_limit_s is not None:
            step_limit = min(step_limit, self.last_run_step_limit_s)
        self.last_run_step_limit_s = step_limit
        rows: list[dict[str, float]] = []
        if tracker is not None and tracker.previous_time_ps is None:
            tracker_summary, cluster_rows = tracker.update(
                self.time_s * 1.0e12,
                self.state.ag_liquid_fraction,
                self.ag_fraction > 0.0,
            )
            rows.append(self._sample_row(0.0, tracker_summary))
        while self.time_s < end_s - 1.0e-24:
            targets = [end_s, next_sample]
            if pending_checkpoints:
                targets.append(pending_checkpoints[0])
            next_target = min(target for target in targets if target > self.time_s + 1.0e-25)
            dt = min(step_limit, next_target - self.time_s)
            if dt <= 1.0e-25:
                continue
            row = self.advance(dt, source_on=source_on, enforce_stability=False)
            at_sample = self.time_s >= next_sample - 1.0e-24
            at_checkpoint = bool(pending_checkpoints) and self.time_s >= pending_checkpoints[0] - 1.0e-24
            at_end = self.time_s >= end_s - 1.0e-24
            if at_sample or at_checkpoint or at_end:
                tracker_summary = None
                if tracker is not None:
                    tracker_summary, cluster_rows = tracker.update(
                        self.time_s * 1.0e12,
                        self.state.ag_liquid_fraction,
                        self.ag_fraction > 0.0,
                    )
                    row.update(tracker_summary)
                row["ag_latent_energy_J"] = self._ag_latent_energy_J()
                rows.append(row)
                if progress:
                    print(
                        f"t={row['time_ps']:.6g} ps  peak={row['max_temperature_K']:.4f} K  "
                        f"Ag liquid={row['ag_liquid_volume_fraction']:.6g}  "
                        f"tracked clusters={row.get('connected_liquid_cluster_count', 0)}",
                        file=sys.stderr,
                        flush=True,
                    )
                if at_sample:
                    next_sample += sample_interval_s
                if at_checkpoint:
                    checkpoint_ps = pending_checkpoints.pop(0) * 1.0e12
                    assert checkpoint_dir is not None
                    checkpoint_name = f"{checkpoint_prefix}_{_time_tag_ps(checkpoint_ps)}.npz"
                    self.save_checkpoint(Path(checkpoint_dir) / checkpoint_name, tracker)
                if at_end:
                    break
        return rows

    def _sample_row(
        self,
        dt_fs: float,
        tracker_summary: dict[str, int | float | list[float] | None] | None = None,
    ) -> dict[str, float]:
        total_energy = float(self.enthalpy_density_J_m3.sum() * self.cell_volume_m3)
        row: dict[str, float] = {
            "time_ps": self.time_s * 1.0e12,
            "dt_fs": dt_fs,
            "max_temperature_K": float(self.state.temperature_K.max()),
            "mean_ag_temperature_K": self._ag_weighted_mean_temperature(),
            "max_ag_liquid_fraction": float(self.state.ag_liquid_fraction.max()),
            "ag_liquid_volume_fraction": self._liquid_ag_volume_fraction(),
            "total_energy_J": total_energy,
            "ag_latent_energy_J": self._ag_latent_energy_J(),
            "injected_energy_J": self.injected_energy_J,
            "energy_balance_error_J": total_energy - self.injected_energy_J,
        }
        if tracker_summary is not None:
            row.update(tracker_summary)
        return row


def _time_tag_ps(time_ps: float) -> str:
    value = f"{time_ps:.9f}".rstrip("0").rstrip(".")
    return value.replace(".", "p") + "ps"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("Cannot write a timeseries with no samples")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_cluster_csv(path: Path, records: list[dict[str, object]]) -> None:
    if not records:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(records[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for record in records:
            row = dict(record)
            for field in ("parent_track_ids", "centroid_A"):
                if field in row:
                    row[field] = json.dumps(row[field], separators=(",", ":"))
            writer.writerow(row)


def longest_sampled_liquid_presence_ps(rows: list[dict[str, float]], sample_interval_ps: float) -> float:
    """Return the longest uninterrupted sampled presence of any liquid cluster."""
    longest = current = 0.0
    previous_time: float | None = None
    tolerance = max(1.0e-12, sample_interval_ps * 1.0e-6)
    for row in rows:
        if row.get("connected_liquid_cluster_count", 0) <= 0:
            current = 0.0
            previous_time = None
            continue
        time_ps = float(row["time_ps"])
        if previous_time is None:
            current = 0.0
        else:
            delta = time_ps - previous_time
            if delta <= sample_interval_ps + tolerance:
                current += delta
            else:
                current = 0.0
        longest = max(longest, current)
        previous_time = time_ps
    return longest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spacing", choices=("5A", "10A"), default="10A")
    parser.add_argument("--geometry", choices=("projected", "native"), default="projected")
    parser.add_argument("--boundary", choices=("A", "B"), required=True)
    parser.add_argument("--heat-ps", type=float, default=1.0)
    parser.add_argument("--cool-ps", type=float, default=1.0)
    parser.add_argument("--sample-ps", type=float, default=0.1)
    parser.add_argument("--tile-xy", type=int, default=1)
    parser.add_argument("--mesh-refinement", type=int, default=1)
    parser.add_argument("--source-integration", choices=("cell-center", "exact-qsgs"), default=None)
    parser.add_argument("--max-dt-fs", type=float, default=None)
    parser.add_argument("--track-liquid-connectivity", action="store_true")
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--liquid-threshold", type=float, default=0.5)
    parser.add_argument("--checkpoint-dir", type=Path, default=None)
    parser.add_argument("--checkpoint-times-ps", type=str, default="")
    parser.add_argument("--checkpoint-prefix", type=str, default="caseA_heat")
    parser.add_argument("--restart-from", type=Path, default=None)
    parser.add_argument("--save-final-checkpoint", action="store_true")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    if args.tile_xy <= 0:
        parser.error("--tile-xy must be a positive integer")
    if args.mesh_refinement <= 0:
        parser.error("--mesh-refinement must be a positive integer")
    if args.geometry != "native" and args.mesh_refinement != 1:
        parser.error("--mesh-refinement applies only to --geometry native")
    if args.max_dt_fs is not None and args.max_dt_fs <= 0.0:
        parser.error("--max-dt-fs must be positive")
    if args.heat_ps < 0.0 or args.cool_ps < 0.0 or args.heat_ps + args.cool_ps <= 0.0:
        parser.error("--heat-ps and --cool-ps must be nonnegative with a positive total duration")
    native_path = DATA_DIR.parent / "qsgs_morphology.npz"
    geometry_path = native_path if args.source_integration == "exact-qsgs" or args.geometry == "native" else None
    if args.track_liquid_connectivity and args.geometry != "native":
        parser.error("--track-liquid-connectivity requires --geometry native")
    if args.geometry == "native" and args.source_integration == "cell-center":
        parser.error("Native geometry requires --source-integration exact-qsgs")
    try:
        checkpoint_times_ps = [
            float(value.strip()) for value in args.checkpoint_times_ps.split(",") if value.strip()
        ]
    except ValueError:
        parser.error("--checkpoint-times-ps must be a comma-separated list of numbers")
    if checkpoint_times_ps and args.checkpoint_dir is None:
        parser.error("--checkpoint-dir is required with --checkpoint-times-ps")
    if args.checkpoint_dir is not None and args.geometry != "native":
        parser.error("Restartable Case A checkpoints require --geometry native")
    if args.restart_from is not None and args.geometry != "native":
        parser.error("--restart-from requires --geometry native")
    projected_path = DATA_DIR / f"qsgs_phase_fractions_{args.spacing}.npz" if args.geometry == "projected" else None
    model = QSGSEnthalpyGrid(
        projected_path,
        args.boundary,
        tile_xy=args.tile_xy,
        source_geometry_npz=geometry_path,
        native_geometry_npz=native_path if args.geometry == "native" else None,
        mesh_refinement=args.mesh_refinement,
    )
    initial_stable_timestep_fs = model.stable_timestep_s() * 1.0e15
    max_dt_s = None if args.max_dt_fs is None else args.max_dt_fs * 1.0e-15
    tracker = (
        LiquidClusterTracker(model.spacing_A, model.size_A, liquid_threshold=args.liquid_threshold)
        if args.track_liquid_connectivity
        else None
    )
    restart_info = None
    if args.restart_from is not None:
        restart_info = model.load_checkpoint(args.restart_from, tracker)
    rows: list[dict[str, object]] = []
    if args.heat_ps > 0.0:
        rows += model.run(
            args.heat_ps,
            source_on=True,
            sample_interval_ps=args.sample_ps,
            maximum_dt_s=max_dt_s,
            track_liquid_connectivity=args.track_liquid_connectivity,
            liquid_threshold=args.liquid_threshold,
            progress=args.progress,
            tracker=tracker,
            checkpoint_times_ps=checkpoint_times_ps,
            checkpoint_dir=args.checkpoint_dir,
            checkpoint_prefix=args.checkpoint_prefix,
        )
    if args.cool_ps > 0.0:
        rows += model.run(
            args.cool_ps,
            source_on=False,
            sample_interval_ps=args.sample_ps,
            maximum_dt_s=max_dt_s,
            track_liquid_connectivity=args.track_liquid_connectivity,
            liquid_threshold=args.liquid_threshold,
            progress=args.progress,
            tracker=tracker,
            checkpoint_times_ps=checkpoint_times_ps,
            checkpoint_dir=args.checkpoint_dir,
            checkpoint_prefix=args.checkpoint_prefix,
        )
    output = args.output or DATA_DIR / f"qsgs_enthalpy_{args.spacing}_BC{args.boundary}.csv"
    write_csv(output, rows)
    mesh_geometry_path = native_path if args.geometry == "native" else projected_path
    assert mesh_geometry_path is not None
    source_geometry_hash = (
        hashlib.sha256(geometry_path.read_bytes()).hexdigest() if geometry_path is not None else None
    )
    summary = {
        "status": "SINGLE_TEMPERATURE_ENTHALPY_EXPLORATORY_RUN",
        "boundary_case": args.boundary,
        "geometry_mode": args.geometry,
        "mesh_geometry_input": str(mesh_geometry_path.relative_to(DATA_DIR.parent)),
        "mesh_geometry_sha256": hashlib.sha256(mesh_geometry_path.read_bytes()).hexdigest(),
        "source_geometry_input": (
            str(geometry_path.relative_to(DATA_DIR.parent)) if geometry_path is not None else None
        ),
        "source_geometry_sha256": source_geometry_hash,
        "spacing_label": args.spacing,
        "tile_xy": args.tile_xy,
        "native_mesh_refinement": args.mesh_refinement if args.geometry == "native" else None,
        "grid_shape_xyz": list(model.temperature_shape),
        "actual_spacing_A": model.spacing_A.tolist(),
        "source_integration_method": model.source_integration_method,
        "stable_timestep_fs_at_initial_state": initial_stable_timestep_fs,
        "maximum_timestep_fs": args.max_dt_fs,
        "liquid_connectivity_threshold": args.liquid_threshold if args.track_liquid_connectivity else None,
        "heat_duration_ps": args.heat_ps,
        "cool_duration_ps": args.cool_ps,
        "start_time_ps": restart_info["time_ps"] if restart_info else 0.0,
        "total_time_ps": model.time_s * 1.0e12,
        "actual_source_power_W": model.actual_source_power_W,
        "injected_energy_J": model.injected_energy_J,
        "final_total_enthalpy_J": float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3),
        "energy_balance_error_J": float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3 - model.injected_energy_J),
        "max_temperature_K": float(model.state.temperature_K.max()),
        "max_temperature_seen_K": float(model.max_temperature_seen_K),
        "mean_ag_temperature_K": model._ag_weighted_mean_temperature(),
        "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
        "max_ag_liquid_fraction": float(model.state.ag_liquid_fraction.max()),
        "ag_normal_boiling_reference_K": AG_NORMAL_BOILING_REFERENCE_K,
        "peak_exceeds_ag_normal_boiling_reference": bool(
            model.max_temperature_seen_K > AG_NORMAL_BOILING_REFERENCE_K
        ),
        "limitations": [
            "Heat is deposited directly into lattice enthalpy; this is not TTM or TTM-MD.",
            "Candidate properties are first-pass values pending source-table archival and citation review.",
            "Boundary A/B are limiting assumptions, not measured Ag/Ti3SiC2 interface heat transfer.",
            "External z faces are insulated and x/y are periodic; no free-surface radiation or convection is included.",
            "The phase model does not include vaporization, so states above the Ag normal-boiling reference are outside its reliable regime.",
            "Adiabatic cooling redistributes energy inside the computational domain and does not represent heat loss to an external sink.",
        ],
        "timeseries_csv": str(output.name),
    }
    summary_path = output.with_suffix(".json")
    if args.track_liquid_connectivity:
        assert tracker is not None
        cluster_path = output.with_name(output.stem + ".clusters.csv")
        write_cluster_csv(cluster_path, tracker.history)
        track_durations: dict[int, float] = {}
        for item in tracker.history:
            track_id = int(item["track_id"])
            track_durations[track_id] = max(
                track_durations.get(track_id, 0.0), float(item["continuous_duration_ps"])
            )
        longest_track_ps = max(track_durations.values(), default=0.0)
        main_track_durations = [
            float(item["continuous_duration_ps"])
            for item in tracker.history
            if bool(item["is_main_cluster"])
        ]
        longest_main_track_ps = max(main_track_durations, default=0.0)
        active_cluster_rows = [
            row for row in rows if row.get("connected_liquid_cluster_count", 0) > 0
        ]
        largest_cluster_row = max(
            active_cluster_rows,
            key=lambda row: row["largest_liquid_cluster_voxels"],
            default=None,
        )
        summary["tracked_liquid_persistence_screen"] = {
            "longest_single_track_sampled_duration_ps": longest_track_ps,
            "longest_track_while_it_was_the_largest_component_ps": longest_main_track_ps,
            "provisional_0p3_ps_gate_met": longest_main_track_ps >= 0.3,
            "provisional_0p5_ps_enhanced_gate_met": longest_main_track_ps >= 0.5,
            "interpretation": (
                "Provisional user-defined screening only. IDs continue through overlap using the documented dominant-overlap rule; "
                "split children get new IDs and merged clusters inherit the largest-overlap parent. This is sampled continuity, "
                "not a universal thermal-stability criterion."
            ),
            "sampling_interval_ps": args.sample_ps,
        }
        summary["liquid_track_onset_times_ps"] = {
            "first_nonzero_liquid_fraction": tracker.first_liquid_time_ps,
            "first_connected_cluster_at_threshold": tracker.first_connected_liquid_time_ps,
            "first_cluster_with_at_least_two_face_connected_voxels": tracker.first_multivoxel_pool_time_ps,
        }
        summary["liquid_tracking"] = {
            "threshold_ag_liquid_fraction": args.liquid_threshold,
            "track_count": len(track_durations),
            "merge_component_event_count": sum(
                "merge" in str(item["event"]) for item in tracker.history
            ),
            "split_parent_event_count": len(
                {
                    (float(item["time_ps"]), int(parent_id))
                    for item in tracker.history
                    if "split" in str(item["event"])
                    for parent_id in item["parent_track_ids"]
                }
            ),
            "cluster_history_csv": cluster_path.name,
            "identity_rule": "Largest-overlap child inherits parent on split; largest-overlap parent supplies ID on merge; all parent IDs remain in the event record.",
        }
        summary["first_detected_connected_liquid_time_ps"] = (
            float(active_cluster_rows[0]["time_ps"]) if active_cluster_rows else None
        )
        summary["connected_liquid_sample_count"] = len(active_cluster_rows)
        summary["max_sampled_connected_liquid_cluster_count"] = max(
            (int(row["connected_liquid_cluster_count"]) for row in rows), default=0
        )
        summary["final_connected_liquid_cluster_count"] = int(
            rows[-1].get("connected_liquid_cluster_count", 0)
        )
        summary["largest_sampled_liquid_cluster"] = (
            {
                "time_ps": float(largest_cluster_row["time_ps"]),
                "voxel_count": int(largest_cluster_row["largest_liquid_cluster_voxels"]),
                "volume_A3": float(largest_cluster_row["largest_liquid_cluster_volume_A3"]),
                "centroid_A": largest_cluster_row["largest_liquid_cluster_centroid_A"],
            }
            if largest_cluster_row is not None
            else None
        )
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    if args.save_final_checkpoint:
        final_checkpoint = output.with_name(output.stem + ".checkpoint.npz")
        model.save_checkpoint(final_checkpoint, tracker)
        summary["final_checkpoint"] = final_checkpoint.name
        summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
