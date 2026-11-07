"""Mass-specific enthalpy and transport functions for the first-pass model."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


DATA_DIR = Path(__file__).resolve().parent
DEFAULT_TABLE = DATA_DIR / "thermophysical_candidates.csv"
DEFAULT_CONFIG = DATA_DIR / "thermophysical_candidates.json"


def _rows_for(table_path: Path, material: str, phase: str) -> list[dict[str, str]]:
    with table_path.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["material"] == material and row["phase"] == phase]
    rows.sort(key=lambda row: float(row["temperature_K"]))
    if not rows:
        raise ValueError(f"No rows for {material}/{phase} in {table_path}")
    return rows


def _column(rows: list[dict[str, str]], name: str) -> np.ndarray:
    return np.asarray([float(row[name]) for row in rows], dtype=float)


def _scalarize(original: np.ndarray | float, value: np.ndarray) -> np.ndarray | float:
    return float(value) if np.asarray(original).ndim == 0 else value


@dataclass(frozen=True)
class AgState:
    temperature_K: np.ndarray | float
    liquid_fraction: np.ndarray | float


class AgEnthalpy:
    """Ag sensible enthalpy plus an isothermal fusion-latent-heat plateau.

    The supplied two solid Cp anchors define a piecewise-linear candidate Cp(T)
    between 300 K and Tm. The first-pass liquid Cp is constant because only one
    liquid anchor is present. These choices are explicit, replaceable inputs,
    not additional measurements.
    """

    def __init__(
        self,
        table_path: Path = DEFAULT_TABLE,
        config_path: Path = DEFAULT_CONFIG,
    ) -> None:
        config = json.loads(config_path.read_text())
        self.table_path = Path(table_path)
        self.molar_mass_kg_mol = float(config["ag"]["molar_mass_kg_mol"])
        self.melting_temperature_K = float(config["ag"]["melting_temperature_K"])
        self.latent_heat_J_kg = float(config["ag"]["latent_heat_fusion_J_kg"])
        self.reference_density_kg_m3 = float(config["ag"]["solid_reference_density_kg_m3"])

        solid = _rows_for(self.table_path, "Ag", "solid")
        liquid = _rows_for(self.table_path, "Ag", "liquid")
        self.solid_temperature_K = _column(solid, "temperature_K")
        self.solid_cp_J_kg_K = _column(solid, "Cp_molar_J_mol_K") / self.molar_mass_kg_mol
        self.solid_k_W_m_K = _column(solid, "thermal_conductivity_W_m_K")
        self.liquid_temperature_K = _column(liquid, "temperature_K")
        self.liquid_cp_J_kg_K = _column(liquid, "Cp_molar_J_mol_K") / self.molar_mass_kg_mol
        self.liquid_k_W_m_K = _column(liquid, "thermal_conductivity_W_m_K")
        self.reference_temperature_K = float(self.solid_temperature_K[0])
        if self.solid_temperature_K[-1] != self.melting_temperature_K:
            raise ValueError("Solid Ag table must end at the configured melting temperature")
        if self.liquid_temperature_K[0] != self.melting_temperature_K:
            raise ValueError("Liquid Ag table must start at the configured melting temperature")
        if self.latent_heat_J_kg <= 0.0 or self.reference_density_kg_m3 <= 0.0:
            raise ValueError("Ag latent heat and reference density must be positive")
        self.solid_h_at_melt_J_kg = float(self._integrate_solid_cp(self.melting_temperature_K))
        self.liquidus_h_J_kg = self.solid_h_at_melt_J_kg + self.latent_heat_J_kg

    def cp_solid_J_kg_K(self, temperature_K: np.ndarray | float) -> np.ndarray | float:
        temperature = np.asarray(temperature_K, dtype=float)
        result = np.interp(temperature, self.solid_temperature_K, self.solid_cp_J_kg_K)
        return _scalarize(temperature_K, result)

    def cp_liquid_J_kg_K(self, temperature_K: np.ndarray | float) -> np.ndarray | float:
        temperature = np.asarray(temperature_K, dtype=float)
        if np.any(temperature < self.melting_temperature_K):
            raise ValueError("Liquid heat capacity is defined only at or above Tm")
        result = np.full_like(temperature, self.liquid_cp_J_kg_K[0], dtype=float)
        return _scalarize(temperature_K, result)

    def conductivity_W_m_K(
        self,
        temperature_K: np.ndarray | float,
        liquid_fraction: np.ndarray | float = 0.0,
    ) -> np.ndarray | float:
        temperature, fraction = np.broadcast_arrays(
            np.asarray(temperature_K, dtype=float), np.asarray(liquid_fraction, dtype=float)
        )
        if np.any((fraction < 0.0) | (fraction > 1.0)):
            raise ValueError("liquid_fraction must lie in [0,1]")
        if np.any(temperature > self.melting_temperature_K) and self.liquid_k_W_m_K.size != 1:
            raise ValueError("Liquid conductivity curve is not yet defined above its single anchor")
        solid_k = np.interp(temperature, self.solid_temperature_K, self.solid_k_W_m_K)
        liquid_k = float(self.liquid_k_W_m_K[0])
        result = np.where(temperature < self.melting_temperature_K, solid_k,
                          (1.0 - fraction) * solid_k + fraction * liquid_k)
        result = np.where(temperature > self.melting_temperature_K, liquid_k, result)
        return float(result) if np.asarray(temperature_K).ndim == 0 and np.asarray(liquid_fraction).ndim == 0 else result

    def cp_equilibrium_J_kg_K(
        self,
        temperature_K: np.ndarray | float,
        liquid_fraction: np.ndarray | float = 0.0,
    ) -> np.ndarray | float:
        temperature, fraction = np.broadcast_arrays(
            np.asarray(temperature_K, dtype=float), np.asarray(liquid_fraction, dtype=float)
        )
        if np.any((fraction < 0.0) | (fraction > 1.0)):
            raise ValueError("liquid_fraction must lie in [0,1]")
        solid_cp = np.interp(temperature, self.solid_temperature_K, self.solid_cp_J_kg_K)
        liquid_cp = float(self.liquid_cp_J_kg_K[0])
        result = np.where(temperature < self.melting_temperature_K, solid_cp,
                          (1.0 - fraction) * solid_cp + fraction * liquid_cp)
        result = np.where(temperature > self.melting_temperature_K, liquid_cp, result)
        return float(result) if np.asarray(temperature_K).ndim == 0 and np.asarray(liquid_fraction).ndim == 0 else result

    def _integrate_solid_cp(self, temperature_K: np.ndarray | float) -> np.ndarray | float:
        temperature = np.asarray(temperature_K, dtype=float)
        if np.any(temperature > self.melting_temperature_K + 1.0e-9):
            raise ValueError("Solid Ag enthalpy is defined only through Tm")
        cumulative = np.zeros_like(self.solid_temperature_K)
        cumulative[1:] = np.cumsum(
            np.diff(self.solid_temperature_K)
            * (self.solid_cp_J_kg_K[:-1] + self.solid_cp_J_kg_K[1:]) / 2.0
        )
        clipped = np.clip(temperature, self.reference_temperature_K, self.melting_temperature_K)
        index = np.clip(np.searchsorted(self.solid_temperature_K, clipped, side="right") - 1,
                        0, len(self.solid_temperature_K) - 2)
        delta = clipped - self.solid_temperature_K[index]
        slope = (self.solid_cp_J_kg_K[index + 1] - self.solid_cp_J_kg_K[index]) / (
            self.solid_temperature_K[index + 1] - self.solid_temperature_K[index]
        )
        integrated = cumulative[index] + self.solid_cp_J_kg_K[index] * delta + 0.5 * slope * delta**2
        below = temperature < self.reference_temperature_K
        integrated = np.where(below,
                              self.solid_cp_J_kg_K[0] * (temperature - self.reference_temperature_K),
                              integrated)
        return _scalarize(temperature_K, integrated)

    def specific_enthalpy_J_kg(
        self,
        temperature_K: np.ndarray | float,
        liquid_fraction: np.ndarray | float = 0.0,
    ) -> np.ndarray | float:
        temperature, fraction = np.broadcast_arrays(
            np.asarray(temperature_K, dtype=float), np.asarray(liquid_fraction, dtype=float)
        )
        if np.any((fraction < 0.0) | (fraction > 1.0)):
            raise ValueError("liquid_fraction must lie in [0,1]")
        tm = self.melting_temperature_K
        hs = np.asarray(self._integrate_solid_cp(np.minimum(temperature, tm)))
        plateau_h = self.solid_h_at_melt_J_kg + fraction * self.latent_heat_J_kg
        liquid_h = self.liquidus_h_J_kg + float(self.liquid_cp_J_kg_K[0]) * (temperature - tm)
        result = np.where(temperature < tm, hs,
                          np.where(temperature > tm, liquid_h, plateau_h))
        if np.any((temperature < tm) & (fraction > 1.0e-12)) or np.any((temperature > tm) & (fraction < 1.0 - 1.0e-12)):
            raise ValueError("Below Tm Ag must be solid; above Tm Ag must be liquid")
        return float(result) if np.asarray(temperature_K).ndim == 0 and np.asarray(liquid_fraction).ndim == 0 else result

    def state_from_specific_enthalpy(self, enthalpy_J_kg: np.ndarray | float) -> AgState:
        enthalpy = np.asarray(enthalpy_J_kg, dtype=float)
        temperature = np.empty_like(enthalpy)
        fraction = np.zeros_like(enthalpy)
        solid = enthalpy < self.solid_h_at_melt_J_kg
        below_reference = solid & (enthalpy < 0.0)
        temperature[below_reference] = self.reference_temperature_K + (
            enthalpy[below_reference] / self.solid_cp_J_kg_K[0]
        )
        needs_inversion = solid & ~below_reference
        if np.any(needs_inversion):
            target = enthalpy[needs_inversion]
            if len(self.solid_temperature_K) == 2:
                cp0 = self.solid_cp_J_kg_K[0]
                slope = (self.solid_cp_J_kg_K[1] - cp0) / (
                    self.solid_temperature_K[1] - self.solid_temperature_K[0]
                )
                if abs(slope) < 1.0e-18:
                    delta = target / cp0
                else:
                    delta = 2.0 * target / (cp0 + np.sqrt(cp0 * cp0 + 2.0 * slope * target))
                temperature[needs_inversion] = self.reference_temperature_K + delta
            else:
                lo = np.full_like(target, self.reference_temperature_K)
                hi = np.full_like(target, self.melting_temperature_K)
                for _ in range(64):
                    mid = 0.5 * (lo + hi)
                    hmid = np.asarray(self._integrate_solid_cp(mid))
                    lower = hmid < target
                    lo = np.where(lower, mid, lo)
                    hi = np.where(lower, hi, mid)
                temperature[needs_inversion] = 0.5 * (lo + hi)

        plateau = (enthalpy >= self.solid_h_at_melt_J_kg) & (enthalpy <= self.liquidus_h_J_kg)
        temperature[plateau] = self.melting_temperature_K
        fraction[plateau] = np.clip(
            (enthalpy[plateau] - self.solid_h_at_melt_J_kg) / self.latent_heat_J_kg,
            0.0,
            1.0,
        )

        liquid = enthalpy > self.liquidus_h_J_kg
        temperature[liquid] = self.melting_temperature_K + (
            enthalpy[liquid] - self.liquidus_h_J_kg
        ) / float(self.liquid_cp_J_kg_K[0])
        fraction[liquid] = 1.0
        if enthalpy.ndim == 0:
            return AgState(float(temperature), float(fraction))
        return AgState(temperature, fraction)


class Ti3SiC2ThermalProperties:
    """Barsoum-based Ti3SiC2 Cp(T), k(T), and analytically integrated enthalpy.

    Cp_molar = 164.4 - 16419/T (T in K), from Barsoum et al. (1999).
    k = 38.6 - 0.0045*T, from the same paper (reported fit R^2=0.86).
    The 2000 K endpoint is an extrapolation, not a measured property value.
    """

    def __init__(
        self,
        table_path: Path = DEFAULT_TABLE,
        config_path: Path = DEFAULT_CONFIG,
    ) -> None:
        config = json.loads(config_path.read_text())
        material = config["ti3sic2"]
        self.molar_mass_kg_mol = float(material["molar_mass_kg_mol"])
        self.density_kg_m3 = float(material["density_kg_m3"])
        self.reference_temperature_K = float(material["reference_temperature_K"])
        self.maximum_temperature_K = float(material["maximum_temperature_K"])
        cp_model = material["cp_model"]
        k_model = material["conductivity_model"]
        self.cp_a_J_mol_K = float(cp_model["a_J_mol_K"])
        self.cp_b_J_mol = float(cp_model["b_J_mol"])
        self.k_intercept_W_m_K = float(k_model["intercept_W_m_K"])
        self.k_slope_W_m_K2 = float(k_model["slope_W_m_K2"])

        rows = _rows_for(Path(table_path), "Ti3SiC2", "solid")
        self.temperature_K = _column(rows, "temperature_K")
        table_cp_molar = _column(rows, "Cp_molar_J_mol_K")
        table_k = _column(rows, "thermal_conductivity_W_m_K")
        if not np.isclose(self.temperature_K[0], self.reference_temperature_K):
            raise ValueError("Ti3SiC2 candidate table must start at the configured reference temperature")
        if not np.isclose(self.temperature_K[-1], self.maximum_temperature_K):
            raise ValueError("Ti3SiC2 candidate table must end at the configured maximum temperature")

        expected_cp = np.asarray(self.cp_molar_J_mol_K(self.temperature_K))
        expected_k = np.asarray(self.conductivity_W_m_K(self.temperature_K))
        if not np.allclose(table_cp_molar, expected_cp, rtol=1.0e-10, atol=1.0e-8):
            raise ValueError("Ti3SiC2 CSV Cp checkpoints do not match the configured Barsoum relation")
        if not np.allclose(table_k, expected_k, rtol=1.0e-10, atol=1.0e-10):
            raise ValueError("Ti3SiC2 CSV conductivity checkpoints do not match the configured Barsoum relation")
        self.cp_mass_J_kg_K = expected_cp / self.molar_mass_kg_mol
        self.conductivity_data_W_m_K = expected_k

    def _check_temperature(self, temperature_K: np.ndarray | float) -> np.ndarray:
        temperature = np.asarray(temperature_K, dtype=float)
        if np.any(temperature <= 0.0):
            raise ValueError("Ti3SiC2 temperature must be positive")
        if np.any(temperature > self.maximum_temperature_K + 1.0e-9):
            raise ValueError(
                f"Ti3SiC2 property model upper limit is {self.maximum_temperature_K:.2f} K"
            )
        return temperature

    def cp_molar_J_mol_K(self, temperature_K: np.ndarray | float) -> np.ndarray | float:
        temperature = self._check_temperature(temperature_K)
        effective = np.maximum(temperature, self.reference_temperature_K)
        result = self.cp_a_J_mol_K - self.cp_b_J_mol / effective
        return _scalarize(temperature_K, result)

    def cp_J_kg_K(self, temperature_K: np.ndarray | float) -> np.ndarray | float:
        result = np.asarray(self.cp_molar_J_mol_K(temperature_K)) / self.molar_mass_kg_mol
        return _scalarize(temperature_K, result)

    def specific_enthalpy_J_kg(self, temperature_K: np.ndarray | float) -> np.ndarray | float:
        """Specific sensible enthalpy relative to 300 K, in J/kg."""
        temperature = self._check_temperature(temperature_K)
        reference = self.reference_temperature_K
        effective = np.maximum(temperature, reference)
        integrated = (
            self.cp_a_J_mol_K * (effective - reference)
            - self.cp_b_J_mol * np.log(effective / reference)
        ) / self.molar_mass_kg_mol
        cp_reference = (
            self.cp_a_J_mol_K - self.cp_b_J_mol / reference
        ) / self.molar_mass_kg_mol
        result = np.where(
            temperature < reference,
            cp_reference * (temperature - reference),
            integrated,
        )
        return _scalarize(temperature_K, result)

    def conductivity_W_m_K(self, temperature_K: np.ndarray | float) -> np.ndarray | float:
        temperature = self._check_temperature(temperature_K)
        result = self.k_intercept_W_m_K + self.k_slope_W_m_K2 * temperature
        if np.any(result <= 0.0):
            raise ValueError("Ti3SiC2 conductivity relation returned a nonpositive value")
        return _scalarize(temperature_K, result)
