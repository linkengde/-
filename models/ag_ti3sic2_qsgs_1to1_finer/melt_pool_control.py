"""Adaptive heat-off controller for externally computed structure metrics.

The controller does not classify liquid atoms. It consumes diagnostic samples
from a validated CNA/MSD/cluster analysis and skeleton-connectivity monitor.
All numerical thresholds below are operational starting values, not material
constants; they must be sensitivity-tested before a production run.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np


class HeatAction(str, Enum):
    CONTINUE = "continue_heating"
    HEAT_OFF = "turn_heat_off_continue_cooling"
    SAFE_STOP = "turn_heat_off_and_stop_run"


@dataclass(frozen=True)
class PoolObservation:
    time_ps: float
    largest_liquid_cluster_atoms: int
    liquid_pool_volume_A3: float
    skeleton_connected: bool
    skeleton_damage_fraction: float
    surface_recession_A: float = 0.0
    escaped_atoms: int = 0


@dataclass(frozen=True)
class HeatDecision:
    action: HeatAction
    reason: str
    stable_duration_ps: float


class StablePoolController:
    """Turn heat off after a persistent, non-growing Ag pool is observed."""

    def __init__(
        self,
        min_liquid_cluster_atoms: int = 500,
        persistence_ps: float = 3.0,
        volume_window_ps: float = 2.0,
        max_relative_volume_drift: float = 0.05,
        max_skeleton_damage_fraction: float = 0.05,
        max_surface_recession_A: float = 3.0,
        max_sample_gap_ps: float = 0.25,
    ) -> None:
        if min_liquid_cluster_atoms < 1:
            raise ValueError("min_liquid_cluster_atoms must be positive")
        if persistence_ps <= 0.0 or volume_window_ps <= 0.0:
            raise ValueError("Persistence and volume windows must be positive")
        if max_relative_volume_drift < 0.0 or max_skeleton_damage_fraction < 0.0:
            raise ValueError("Stability tolerances cannot be negative")
        self.min_cluster = int(min_liquid_cluster_atoms)
        self.persistence_ps = float(persistence_ps)
        self.volume_window_ps = float(volume_window_ps)
        self.max_volume_drift = float(max_relative_volume_drift)
        self.max_skeleton_damage = float(max_skeleton_damage_fraction)
        self.max_surface_recession = float(max_surface_recession_A)
        self.max_sample_gap = float(max_sample_gap_ps)
        self._samples: list[PoolObservation] = []
        self._candidate_start: float | None = None
        self._last_time: float | None = None
        self.heat_is_on = True

    def _volume_plateau(self, time_ps: float) -> bool:
        selected = [
            item for item in self._samples
            if time_ps - self.volume_window_ps - 1.0e-12 <= item.time_ps <= time_ps + 1.0e-12
        ]
        if len(selected) < 2 or selected[-1].time_ps - selected[0].time_ps < self.volume_window_ps * 0.95:
            return False
        times = np.asarray([item.time_ps for item in selected], dtype=float)
        volumes = np.asarray([item.liquid_pool_volume_A3 for item in selected], dtype=float)
        mean_volume = float(np.mean(volumes))
        if mean_volume <= 0.0:
            return False
        slope = float(np.polyfit(times - times[0], volumes, 1)[0])
        drift = abs(slope) * self.volume_window_ps / mean_volume
        return drift <= self.max_volume_drift

    def update(self, observation: PoolObservation) -> HeatDecision:
        if not self.heat_is_on:
            return HeatDecision(HeatAction.HEAT_OFF, "heat was already turned off", 0.0)
        if not np.isfinite(observation.time_ps) or observation.time_ps < 0.0:
            raise ValueError("Observation time must be finite and non-negative")
        if self._last_time is not None and observation.time_ps <= self._last_time:
            raise ValueError("Pool observations must have strictly increasing times")
        if not np.isfinite(observation.liquid_pool_volume_A3) or observation.liquid_pool_volume_A3 < 0.0:
            raise ValueError("Pool volume must be finite and non-negative")
        if observation.largest_liquid_cluster_atoms < 0 or observation.escaped_atoms < 0:
            raise ValueError("Atom counts cannot be negative")
        if not np.isfinite(observation.skeleton_damage_fraction) or not 0.0 <= observation.skeleton_damage_fraction <= 1.0:
            raise ValueError("Skeleton damage fraction must be in [0, 1]")
        if not np.isfinite(observation.surface_recession_A):
            raise ValueError("Surface recession must be finite")

        self._last_time = observation.time_ps
        self._samples.append(observation)
        cutoff = observation.time_ps - max(self.persistence_ps, self.volume_window_ps) - self.max_sample_gap
        self._samples = [item for item in self._samples if item.time_ps >= cutoff]

        safety_reasons = []
        if observation.escaped_atoms > 0:
            safety_reasons.append(f"{observation.escaped_atoms} atom(s) escaped")
        if not observation.skeleton_connected:
            safety_reasons.append("Ti3SiC2 skeleton lost connectivity")
        if observation.skeleton_damage_fraction > self.max_skeleton_damage:
            safety_reasons.append("skeleton damage exceeded limit")
        if observation.surface_recession_A > self.max_surface_recession:
            safety_reasons.append("surface recession exceeded limit")
        if safety_reasons:
            self.heat_is_on = False
            return HeatDecision(
                HeatAction.SAFE_STOP,
                "; ".join(safety_reasons),
                0.0,
            )

        pool_present = (
            observation.largest_liquid_cluster_atoms >= self.min_cluster
            and observation.liquid_pool_volume_A3 > 0.0
        )
        continuous = True
        if len(self._samples) >= 2:
            gaps = np.diff([item.time_ps for item in self._samples])
            continuous = bool(np.all(gaps <= self.max_sample_gap))
        if not pool_present or not continuous:
            self._candidate_start = observation.time_ps if pool_present else None
            return HeatDecision(HeatAction.CONTINUE, "pool persistence/size gate not met", 0.0)

        if self._candidate_start is None:
            self._candidate_start = observation.time_ps
        duration = observation.time_ps - self._candidate_start
        if duration >= self.persistence_ps and self._volume_plateau(observation.time_ps):
            self.heat_is_on = False
            return HeatDecision(
                HeatAction.HEAT_OFF,
                "connected liquid Ag cluster persisted and pool-volume drift is within tolerance",
                duration,
            )
        return HeatDecision(HeatAction.CONTINUE, "stable-pool window is still accumulating", duration)
