"""Small deterministic checks for the isolated finite-interface kernel."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


PHASE_FIELD_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PHASE_FIELD_DIR))

from finite_interface_prototype import FiniteInterfaceFV  # noqa: E402
from thermophysical_model import AgEnthalpy, Ti3SiC2ThermalProperties  # noqa: E402


class FiniteInterfaceConductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        ag = AgEnthalpy()
        sc = Ti3SiC2ThermalProperties()
        cls.k_ag = float(ag.conductivity_W_m_K(300.0))
        cls.k_sc = float(sc.conductivity_W_m_K(300.0))
        cls.capacity_ag = ag.reference_density_kg_m3 * float(ag.cp_solid_J_kg_K(300.0))
        cls.capacity_sc = sc.density_kg_m3 * float(sc.cp_J_kg_K(300.0))

    def _bilayer(
        self,
        interface_conductance: float,
        n_ag: int = 4,
        n_sc: int = 4,
        dx: float = 1.0e-6,
    ) -> tuple[FiniteInterfaceFV, float, float]:
        area = 1.0e-12
        spacing = (dx, 1.0e-6, 1.0e-6)
        phase = np.zeros((n_ag + n_sc, 1, 1), dtype=np.int8)
        phase[n_ag:, :, :] = 1
        conductivity = np.where(phase == 0, self.k_ag, self.k_sc)
        capacity = np.where(phase == 0, self.capacity_ag, self.capacity_sc)
        model = FiniteInterfaceFV(
            phase,
            spacing,
            conductivity,
            capacity,
            interface_conductance,
            dirichlet_boundaries={(0, "low"): 400.0, (0, "high"): 300.0},
        )
        return model, dx, area

    def test_interface_face_uses_series_bulk_and_kapitza_resistances(self) -> None:
        interface_g = 2.0e8  # benchmark fixture only; not a material estimate
        model, dx, area = self._bilayer(interface_g)
        expected = area / (dx / (2.0 * self.k_ag) + 1.0 / interface_g + dx / (2.0 * self.k_sc))
        np.testing.assert_allclose(model.interface_conductances_W_K, [expected], rtol=2.0e-15, atol=0.0)

    def test_steady_bilayer_matches_series_resistance_for_zero_finite_and_perfect_contact(self) -> None:
        area = 1.0e-12
        delta_t = 100.0
        physical_phase_thickness = 4.0e-6
        for cells_per_phase in (4, 8, 16):
            dx = physical_phase_thickness / cells_per_phase
            for interface_g in (0.0, 2.0e8, float("inf")):
                with self.subTest(cells_per_phase=cells_per_phase, interface_conductance=interface_g):
                    model, _, _ = self._bilayer(interface_g, cells_per_phase, cells_per_phase, dx)
                    temperature = model.steady_state()
                    left_input = 2.0 * self.k_ag * area / dx * (400.0 - temperature[0, 0, 0])
                    right_output = 2.0 * self.k_sc * area / dx * (temperature[-1, 0, 0] - 300.0)
                    if interface_g == 0.0:
                        expected_power = 0.0
                        self.assertAlmostEqual(float(np.ptp(temperature[:cells_per_phase, 0, 0])), 0.0, places=10)
                        self.assertAlmostEqual(float(np.ptp(temperature[cells_per_phase:, 0, 0])), 0.0, places=10)
                    else:
                        interface_resistance = 0.0 if np.isinf(interface_g) else 1.0 / interface_g
                        total_areal_resistance = (
                            physical_phase_thickness / self.k_ag
                            + interface_resistance
                            + physical_phase_thickness / self.k_sc
                        )
                        expected_power = delta_t * area / total_areal_resistance
                    self.assertAlmostEqual(left_input, expected_power, delta=max(abs(expected_power) * 2e-11, 1e-15))
                    self.assertAlmostEqual(right_output, expected_power, delta=max(abs(expected_power) * 2e-11, 1e-15))
                    residual = model.net_power_W(temperature)
                    self.assertLess(float(np.max(np.abs(residual))), max(abs(expected_power) * 2e-11, 1e-14))

    def test_insulated_transient_conserves_total_energy_with_candidate_phase_capacities(self) -> None:
        phase = np.array([0, 0, 1, 1], dtype=np.int8).reshape((4, 1, 1))
        conductivity = np.where(phase == 0, self.k_ag, self.k_sc)
        capacity = np.where(phase == 0, self.capacity_ag, self.capacity_sc)
        model = FiniteInterfaceFV(
            phase,
            (1.0e-6, 1.0e-6, 1.0e-6),
            conductivity,
            capacity,
            2.0e8,  # synthetic benchmark value; no physical interface fit is implied
        )
        initial = np.array([400.0, 400.0, 300.0, 300.0]).reshape(phase.shape)
        dt = 0.2 * model.stable_timestep_s()
        advanced, ledger = model.step(initial, dt)
        self.assertLess(abs(ledger.residual_J), 2.0e-14 * ledger.energy_before_J)
        self.assertAlmostEqual(model.energy_J(advanced), model.energy_J(initial), delta=2.0e-14 * ledger.energy_before_J)
        self.assertLess(float(advanced[1, 0, 0]), 400.0)
        self.assertGreater(float(advanced[2, 0, 0]), 300.0)
        self.assertAlmostEqual(float(np.sum(model.net_power_W(initial))), 0.0, delta=1e-18)

    def test_zero_interface_conductance_blocks_phase_to_phase_exchange(self) -> None:
        phase = np.array([0, 0, 1, 1], dtype=np.int8).reshape((4, 1, 1))
        conductivity = np.where(phase == 0, self.k_ag, self.k_sc)
        capacity = np.where(phase == 0, self.capacity_ag, self.capacity_sc)
        model = FiniteInterfaceFV(
            phase,
            (1.0e-6, 1.0e-6, 1.0e-6),
            conductivity,
            capacity,
            0.0,
        )
        initial = np.array([400.0, 400.0, 300.0, 300.0]).reshape(phase.shape)
        advanced, _ = model.step(initial, 0.25 * model.stable_timestep_s())
        np.testing.assert_array_equal(advanced, initial)
        np.testing.assert_array_equal(model.interface_conductances_W_K, [0.0])

    def test_explicit_step_rejects_timestep_above_stability_bound(self) -> None:
        model, _, _ = self._bilayer(2.0e8)
        with self.assertRaisesRegex(ValueError, "stability bound"):
            model.step(np.full(model.shape, 350.0), 2.0 * model.stable_timestep_s())


if __name__ == "__main__":
    unittest.main()
