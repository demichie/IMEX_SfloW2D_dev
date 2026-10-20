"""Independent identities for the current KT reference, including negative controls."""
from dataclasses import replace
import unittest

import numpy as np

from compare_fields import kt


class KTReferenceTests(unittest.TestCase):
    def test_lake_at_rest_is_preserved(self):
        y, x = np.mgrid[:9, :12]
        vertices = 0.1*x+0.03*y+0.02*np.sin(x*y)
        h = 3-kt.VertexBed(vertices).B_center
        state = np.zeros(h.shape+(4,))
        state[..., 0] = 1000*h
        np.testing.assert_array_equal(kt.rhs(state, vertices, .4, .4, kt.Parameters()), 0)

    def test_path_integral_matches_constant_coefficient_identity(self):
        p = kt.Parameters()
        hm, hp = np.array([1., .2]), np.array([2., .7])
        em, ep = np.array([4., 3.]), np.array([4.5, 2.])
        gamma, gravity = np.full(2, 9800.), np.full(2, .8)
        actual = kt.hydrostatic_path(hm, hp, gamma, gamma, em, ep, gravity, gravity, p)
        np.testing.assert_allclose(actual, -gravity*gamma*.5*(hm+hp)*(ep-em), rtol=4e-16)

    def test_open_boundary_budget_is_not_a_closed_domain_assumption(self):
        state = np.zeros((5, 8, 4))
        state[..., 0], state[..., 1] = 1000., 250.
        vertices = np.zeros((6, 9))
        dt, dx, dy = .001, .5, .75
        final, _, budget = kt.advance(state, vertices, dx, dy, [dt], kt.Parameters(),
                                      return_diagnostics=True)
        expected = -dt*1000*.25*5*dy
        self.assertAlmostEqual(budget['boundary_mass_transfer'], expected, places=13)
        change = (final[..., 0].sum()-state[..., 0].sum())*dx*dy
        self.assertAlmostEqual(change, expected, places=11)
        # Omitting the outgoing flux would incorrectly reject a conservative update.
        self.assertGreater(abs(change), 1e-3)

    def test_axis_transposition_and_collapsed_dimension(self):
        x = np.arange(13)
        vertices = np.stack((.1*np.sin(x), .1*np.sin(x)))
        h = 1+.1*np.cos(np.arange(12))
        state = np.zeros((1, 12, 4))
        state[0, :, 0], state[0, :, 1] = 1000*h, 200*h
        transposed = np.swapaxes(state, 0, 1).copy()
        transposed[..., [1, 2]] = transposed[..., [2, 1]]
        for mode in ('h_x', 'eta_x', 'minmod', 'froude_blend'):
            p = replace(kt.Parameters(), reconstruction=mode)
            original = kt.rhs(state, vertices, .4, .4, p)
            rotated = kt.rhs(transposed, vertices.T, .4, .4, p)
            rotated = np.swapaxes(rotated, 0, 1)
            rotated[..., [1, 2]] = rotated[..., [2, 1]]
            np.testing.assert_allclose(original, rotated, rtol=2e-15, atol=1e-12)

    def test_uniform_solid_mass_fraction_is_passively_transported(self):
        y, x = np.mgrid[:6, :9]
        h, fraction = 1+.1*np.cos(x+y), .25
        rho = 1/((1-fraction)/1000+fraction/2500)
        state = np.zeros(h.shape+(4,))
        state[..., 0] = rho*h
        state[..., 1] = state[..., 0]*(.1+.02*np.sin(x))
        state[..., 2] = state[..., 0]*.03*np.cos(y)
        state[..., 3] = fraction*state[..., 0]
        residual = kt.rhs(state, np.zeros((7, 10)), .5, .5, kt.Parameters())
        np.testing.assert_allclose(residual[..., 3], fraction*residual[..., 0],
                                   rtol=2e-14, atol=2e-12)

    def test_reference_rejects_invalid_raw_states(self):
        state = np.zeros((3, 4, 4))
        state[..., 0] = 1000
        for value in (-1., np.nan):
            bad = state.copy()
            bad[1, 1, 0] = value
            with self.assertRaises(ValueError):
                kt.rhs(bad, np.zeros((4, 5)), 1., 1., kt.Parameters())
        bad = state.copy()
        bad[1, 1, 3] = 1001
        with self.assertRaises(ValueError):
            kt.rhs(bad, np.zeros((4, 5)), 1., 1., kt.Parameters())


if __name__ == '__main__':
    unittest.main()
