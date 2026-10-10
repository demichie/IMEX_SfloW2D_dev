"""Negative controls for geometry/rotation and the independent 2D oracle."""
import unittest
from dataclasses import replace

import numpy as np

import check_cases as c


class CheckerTests(unittest.TestCase):
    def test_rotation_vector_and_inverse(self):
        q = np.zeros((8, 8, 5)); q[2, 3] = (1000, 200, 300, 1e8, 0)
        rotated = c.rotate(q)
        np.testing.assert_array_equal(rotated[4, 2], (1000, 300, -200, 1e8, 0))
        np.testing.assert_array_equal(c.rotate(q, 4), q)
        with self.assertRaises(AssertionError):
            c.equilibrium.field_error("wrong vector rotation", np.rot90(q), rotated, q, c.ROUND)

    def test_isotropy_detects_square_despite_exact_D4(self):
        q = c.equilibrium.liquid_state(np.zeros((48, 48)))
        q[16:32, 16:32] = c.equilibrium.liquid_state(np.ones((16, 16)))
        np.testing.assert_array_equal(q, c.rotate(q))
        self.assertGreater(c.circular_anisotropy(q, 12), c.CONTRACT["circular"]["maximum_fourfold_anisotropy"])
        with self.assertRaises(AssertionError):
            c.circular_anisotropy(np.zeros_like(q), 12)

    def test_flat_reference_is_stationary_all_tableaux(self):
        initial = c.equilibrium.liquid_state(np.ones((8, 8)))
        for n in c.CONTRACT["stages"]:
            result = c.reference(initial, np.zeros((9, 9)), 1, 1, [0.001]*2, n)
            np.testing.assert_array_equal(result["q"], initial)
            self.assertEqual(result["minimum_raw_h"], 1)
        bad = initial.copy(); bad[4, 4, 0] = -1
        with self.assertRaises(AssertionError):
            c.reference(bad, np.zeros((9, 9)), 1, 1, [0.001], 2)

    def test_seeded_reuse_and_common_collar(self):
        for name in c.CONTRACT["seeded"]["fixtures"]:
            vertices, initial, _, _ = c.seeded(name)
            source, _, q = c.spatial.fixture(name)
            pad = c.CONTRACT["seeded"]["pad_cells"]
            np.testing.assert_array_equal(vertices[pad:-pad, pad:-pad], source)
            np.testing.assert_array_equal(c.compact(initial)[pad:-pad, pad:-pad], q)
            self.assertTrue(np.all(initial[0, :, 1:3] == 0))
            self.assertTrue(np.all(initial[-1, :, 1:3] == 0))

    def test_positive_cell_mass_below_face_cutoff_is_not_discarded(self):
        initial = c.equilibrium.liquid_state(np.zeros((8, 8)))
        initial[4, 4] = c.equilibrium.liquid_state(np.array(5e-11))
        result = c.reference(initial, np.zeros((9, 9)), 1, 1, [0.001], 2)
        self.assertGreater(float(result["q"][..., 0].sum()), 0)
        self.assertAlmostEqual(float(result["q"][..., 0].sum()),
                               float(initial[..., 0].sum()), places=15)

    def test_final_dry_projection_retains_mass_and_raw_momentum_evidence(self):
        initial = c.equilibrium.liquid_state(np.full((8, 8), 2.2e-15))
        initial[..., 1] = 5e-6
        initial[..., 2] = 8e-6
        result = c.reference(initial, np.zeros((9, 9)), 1, 1, [0.001], 2)
        np.testing.assert_array_equal(result["q"][..., 0], initial[..., 0])
        np.testing.assert_array_equal(result["q"][..., 3:], result["raw"][..., 3:])
        np.testing.assert_array_equal(result["q"][..., 1:3], 0)
        self.assertGreater(float(np.max(np.abs(result["raw"][..., 1:3]))), 0)
        self.assertGreater(result["repairs"][0], 0)
        np.testing.assert_array_equal(result["known"][0], c.lift(c.compact(initial)))

    def test_stationary_guard_does_not_remove_nonzero_motion(self):
        initial = c.equilibrium.liquid_state(np.ones((8, 8)))
        params = c.core.SolverParams(rest_flux_tol=0)
        bed = c.core.continuous_q1_bed_geometry(np.zeros((9, 9)))
        cache = c.core.prepare_stage_cache2d(c.compact(initial), bed, 0, params,
                                           q_is_sanitized=True, dx=1, dy=1)
        c.final_dry_face_momenta(cache, params, c.core)
        rhs = np.zeros((8, 8, 4)); rhs[..., 1:3] = 1e-15
        c.stationary_roundoff_policy(rhs, cache, 1, 1, params, c.core)
        np.testing.assert_array_equal(rhs, 0)
        initial[..., 1] = 1e-30
        cache = c.core.prepare_stage_cache2d(c.compact(initial), bed, 0, params,
                                           q_is_sanitized=True, dx=1, dy=1)
        c.final_dry_face_momenta(cache, params, c.core)
        rhs[..., 1:3] = 1e-15; original = rhs.copy()
        c.stationary_roundoff_policy(rhs, cache, 1, 1, params, c.core)
        np.testing.assert_array_equal(rhs, original)

    def test_external_ghost_copies_final_trace_not_cell_average(self):
        vertices, initial, dx, dy = c.ramp(45)
        params = c.core.SolverParams(rest_flux_tol=0)
        cache = c.core.prepare_stage_cache2d(c.compact(initial),
            c.core.continuous_q1_bed_geometry(vertices), 0, params,
            q_is_sanitized=True, dx=dx, dy=dy)
        # Inject an unmistakable old ghost only in this negative control.
        data = cache['x'][0]['cache1d']['data']
        data['qL_face'][0] = (1, 2, 0)
        c.final_dry_face_momenta(cache, params, c.core)
        for axis in ('x', 'y'):
            for entry in cache[axis]:
                data = entry['cache1d']['data']
                for prefix in ('q', 'u', 'hu', 'eta'):
                    np.testing.assert_array_equal(data[f'{prefix}L_face'][0], data[f'{prefix}R_face'][0])
                    np.testing.assert_array_equal(data[f'{prefix}R_face'][-1], data[f'{prefix}L_face'][-1])

    def test_all_LS_derivatives_share_tensor_boundary_strip(self):
        vertices, initial, dx, dy = c.mixed()
        bed = c.core.continuous_q1_bed_geometry(vertices)
        geometry = c.production_slope_geometry(bed, dx, dy, c.core)
        original = c.core.build_slope_geometry_2d(bed, dx, dy)
        for name in ('Bx_center', 'By_center', 'Bxx_center', 'Bxy_center', 'Byy_center'):
            field = getattr(geometry, name)
            np.testing.assert_array_equal(field[2:-2, 2:-2], getattr(original, name)[2:-2, 2:-2])
            np.testing.assert_array_equal(field[0], field[2])
            np.testing.assert_array_equal(field[:, 0], field[:, 2])

    def test_dry_auxiliary_velocity_cannot_change_wet_tangential_candidates(self):
        params = c.core.SolverParams(eps_sing=1e-8)
        row = np.zeros((8, 4)); row[:4, 0] = 1000; row[:4, 2] = 2
        expected = c.auxiliary_tangential_faces(row, 'x', params, c.core)
        perturbed = row.copy(); perturbed[4, 0] = 1.7e-14
        # A negligible positive mass can carry a sizable desingularized
        # auxiliary velocity. Neither raw cell mass nor momentum is removed.
        perturbed[4, 2] = 1e-5
        original = perturbed.copy()
        observed = c.auxiliary_tangential_faces(perturbed, 'x', params, c.core)
        for actual, reference in zip(observed, expected):
            np.testing.assert_array_equal(actual, reference)
        np.testing.assert_array_equal(perturbed, original)
        self.assertGreater(abs(c.core._tangential_velocity_slice(perturbed, 'x', params)[4]), 0)
        old_faces = c.core._reconstruct_scalar_faces(
            c.core._tangential_velocity_slice(perturbed, 'x', params), params)
        self.assertTrue(any(np.any(old != new) for old, new in zip(old_faces, observed)))

    def test_curvature_guard_keeps_wet_source_and_dry_conservative_fields(self):
        vertices, initial, dx, dy = c.mixed()
        params = c.core.SolverParams(curvature_term=True, slope_correction=True,
                                   rest_flux_tol=0)
        bed = c.core.continuous_q1_bed_geometry(vertices)
        q = c.compact(initial)
        q[20, 20] = (params.rho_c*params.dry_h_tol/2, 1e-5, -2e-5, 0)
        original = q.copy()
        cache = c.core.prepare_stage_cache2d(q, bed, 0, params,
            q_is_sanitized=True, dx=dx, dy=dy,
            slope_geometry=c.production_slope_geometry(bed, dx, dy, c.core))
        base = c.core.rhs2d(q, bed, dx, dy, 0,
            replace(params, curvature_term=False), stage_cache=cache)
        corrected = c.production_curvature_rhs(q, bed, dx, dy, 0, params, cache, c.core)
        historical = c.core.rhs2d(q, bed, dx, dy, 0, params, stage_cache=cache)
        wet = c.core.primitive2d(q, params)["h"] > params.dry_h_tol
        np.testing.assert_array_equal(corrected[~wet], base[~wet])
        np.testing.assert_array_equal(corrected[wet], historical[wet])
        self.assertGreater(float(np.max(np.abs(historical[~wet]-base[~wet]))), 0)
        np.testing.assert_array_equal(q, original)


if __name__ == "__main__":
    unittest.main()
