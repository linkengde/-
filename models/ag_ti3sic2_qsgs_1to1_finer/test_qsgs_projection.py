import sys
import unittest
from pathlib import Path

import numpy as np

PHASE_DIR = Path(__file__).resolve().parent / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_DIR))

from qsgs_projection import _overlap_weights, project_phase_mask  # noqa: E402
from analyze_projected_connectivity import _component_metrics  # noqa: E402
from audit_coarse_ag_face_connectivity import audit_face_connectivity  # noqa: E402
from qsgs_face_contacts import exact_face_contact_fractions  # noqa: E402
from generate_nested_qsgs_grids import coarsen_face_contacts_2x  # noqa: E402


class QSGSProjectionTests(unittest.TestCase):
    def test_overlap_weights_cover_every_target_cell(self):
        weights = _overlap_weights(7, 4, 3.25)
        np.testing.assert_allclose(weights.sum(axis=1), 1.0, atol=1.0e-12, rtol=0.0)
        self.assertTrue(np.all(weights >= 0.0))

    def test_uniform_projection_is_exact(self):
        skeleton = np.zeros((4, 3, 2), dtype=bool)
        result = project_phase_mask(skeleton, np.array([8.0, 6.0, 4.0]), 5.0)
        np.testing.assert_allclose(result.ag_fraction, 1.0)
        np.testing.assert_allclose(result.skeleton_fraction, 0.0)

    def test_half_volume_is_conserved_in_a_single_target_cell(self):
        skeleton = np.zeros((2, 2, 2), dtype=bool)
        skeleton[0] = True
        result = project_phase_mask(skeleton, np.array([2.0, 2.0, 2.0]), 3.0)
        self.assertEqual(result.ag_fraction.shape, (1, 1, 1))
        self.assertAlmostEqual(float(result.ag_fraction.item()), 0.5)
        self.assertAlmostEqual(float(result.skeleton_fraction.item()), 0.5)

    def test_repository_qsgs_volume_is_conserved_at_both_target_resolutions(self):
        model_dir = Path(__file__).resolve().parent
        with np.load(model_dir / "qsgs_morphology.npz", allow_pickle=False) as saved:
            skeleton = saved["skeleton"]
            size_A = saved["size_A"]
        source_fraction = float(skeleton.mean())
        for requested_spacing_A in (5.0, 10.0):
            with self.subTest(spacing=requested_spacing_A):
                result = project_phase_mask(skeleton, size_A, requested_spacing_A)
                self.assertAlmostEqual(float(result.skeleton_fraction.mean()), source_fraction, places=12)
                self.assertAlmostEqual(float(result.ag_fraction.mean()), 1.0 - source_fraction, places=12)
                np.testing.assert_allclose(
                    result.ag_fraction + result.skeleton_fraction,
                    1.0,
                    atol=1.0e-12,
                    rtol=0.0,
                )

    def test_connectivity_audit_wraps_xy_but_keeps_z_open(self):
        periodic_pair = np.zeros((4, 3, 3), dtype=bool)
        periodic_pair[0, 1, 1] = True
        periodic_pair[-1, 1, 1] = True
        periodic_metrics = _component_metrics(periodic_pair)
        self.assertEqual(periodic_metrics["component_count"], 1)
        self.assertEqual(periodic_metrics["largest_component_fraction"], 1.0)

        open_z_pair = np.zeros((4, 3, 3), dtype=bool)
        open_z_pair[1, 1, 0] = True
        open_z_pair[1, 1, -1] = True
        open_z_metrics = _component_metrics(open_z_pair)
        self.assertEqual(open_z_metrics["component_count"], 2)

    def test_coarse_face_audit_detects_link_without_source_ag_contact(self):
        source_ag = np.zeros((6, 1, 1), dtype=bool)
        source_ag[[1, 4], 0, 0] = True
        coarse_ag = np.array([[[1.0 / 3.0]], [[1.0 / 3.0]]])
        metrics = audit_face_connectivity(
            source_ag,
            coarse_ag,
            np.array([6.0, 1.0, 1.0]),
            np.array([3.0, 1.0, 1.0]),
        )
        x_result = metrics["axis_results"][0]
        self.assertEqual(x_result["coarse_positive_conductance_faces"], 2)
        self.assertEqual(x_result["faces_with_actual_qsgs_ag_contact"], 0)
        self.assertEqual(x_result["coarse_positive_faces_without_actual_ag_contact"], 2)
        self.assertEqual(metrics["coarse_solver_connectivity"]["component_count"], 1)
        self.assertEqual(metrics["direct_source_ag_contact_connectivity"]["component_count"], 2)

    def test_exact_face_contact_fraction_blocks_coarse_false_link(self):
        source_ag = np.zeros((6, 1, 1), dtype=bool)
        source_ag[[1, 4], 0, 0] = True
        coarse_ag = np.array([[[1.0 / 3.0]], [[1.0 / 3.0]]])
        contact_x, _, _ = exact_face_contact_fractions(
            source_ag,
            coarse_ag,
            np.array([6.0, 1.0, 1.0]),
            np.array([3.0, 1.0, 1.0]),
        )
        np.testing.assert_array_equal(contact_x, np.zeros_like(coarse_ag))

    def test_exact_face_contact_fraction_preserves_continuous_ag_link(self):
        source_ag = np.ones((4, 2, 2), dtype=bool)
        coarse_ag = np.ones((2, 1, 1), dtype=float)
        contact_x, contact_y, _ = exact_face_contact_fractions(
            source_ag,
            coarse_ag,
            np.array([4.0, 2.0, 2.0]),
            np.array([2.0, 2.0, 2.0]),
        )
        np.testing.assert_allclose(contact_x, 1.0)
        np.testing.assert_allclose(contact_y, 1.0)

    def test_nested_phase_projection_coarsens_exactly_to_coarse_projection(self):
        rng = np.random.default_rng(20261026)
        skeleton = rng.random((9, 7, 5)) < 0.43
        size_A = np.asarray([30.0, 31.0, 16.0])
        coarse = project_phase_mask(skeleton, size_A, 10.0, (3, 4, 2))
        fine = project_phase_mask(skeleton, size_A, 5.0, (6, 8, 4))
        fine_ag_coarsened = fine.ag_fraction.reshape(3, 2, 4, 2, 2, 2).mean(axis=(1, 3, 5))
        np.testing.assert_allclose(fine_ag_coarsened, coarse.ag_fraction, rtol=0.0, atol=2.0e-15)

    def test_nested_face_contact_areas_coarsen_exactly(self):
        rng = np.random.default_rng(20261027)
        source_ag = rng.random((8, 8, 8)) < 0.6
        size_A = np.asarray([8.0, 8.0, 8.0])
        coarse = project_phase_mask(~source_ag, size_A, 4.0, (2, 2, 2))
        fine = project_phase_mask(~source_ag, size_A, 2.0, (4, 4, 4))
        coarse_contacts = exact_face_contact_fractions(
            source_ag, coarse.ag_fraction, size_A, coarse.spacing_A
        )
        fine_contacts = exact_face_contact_fractions(
            source_ag, fine.ag_fraction, size_A, fine.spacing_A
        )
        for actual, expected in zip(
            coarsen_face_contacts_2x(fine_contacts), coarse_contacts
        ):
            np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2.0e-15)


if __name__ == "__main__":
    unittest.main()
