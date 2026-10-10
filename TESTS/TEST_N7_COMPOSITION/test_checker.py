"""Independent construction and corruption controls for the N7-B comparator."""
import unittest
import hashlib
import json
import numpy as np
import check_cases as checker


class CheckerTests(unittest.TestCase):
    def test_stationary_guard_contract_and_reference_pin(self):
        guard=json.loads((checker.HERE/'roundoff_guard_contract.json').read_text())
        reference=checker.HERE.parent/'TEST_HYDROSTATIC_ROUNDOFF/reference/hyperbolic_2d.f90'
        self.assertEqual(hashlib.sha256(reference.read_bytes()).hexdigest(),guard['baseline_module_sha256'])
        self.assertEqual(guard['correction']['roundoff_multiplier'],64)
        self.assertEqual(guard['controls']['expected_cases'],36)
        self.assertGreaterEqual(guard['controls']['transverse_cells'],5)
        self.assertEqual(guard['controls']['velocities'],[1e-20,1e-12])

    def test_density_nullspace_and_pressure_equilibrium(self):
        profile = np.linspace(0, 1, 24).reshape(3, 8)
        for name in checker.CONTRACT["closures"]:
            ys, base, _ = checker.composition(name, profile)
            rho, _, remainder = checker.properties(name, ys, 300)
            self.assertLess(np.ptp(rho)/rho.max(), checker.ROUND)
            self.assertGreater(ys.min(), 0)
            self.assertGreater(remainder.min(), 0)
            ys, base, _ = checker.composition(name, profile, neutral=False)
            rho, _, _ = checker.properties(name, ys, 300)
            rho0, _, _ = checker.properties(name, base, 300)
            h = np.sqrt((rho0-checker.RHO_A)/(rho-checker.RHO_A))
            pressure = (rho-checker.RHO_A)*h*h
            self.assertLess(np.ptp(pressure)/pressure.max(), checker.ROUND)
            self.assertGreater(np.ptp(h), 1e-5)
            # Holding eta/h fixed with variable density does NOT satisfy equilibrium.
            self.assertGreater(np.ptp(rho-checker.RHO_A), 1e-3)

    def test_raw_thermal_and_fraction_decode(self):
        for name in checker.CONTRACT["closures"]:
            ys, _, _ = checker.composition(name, np.ones((1, 8))*0.5)
            state = checker.conservative(name, ys, np.ones((1, 8)))
            fractions, temperature, _, h, remainder = checker.decode(name, state)
            np.testing.assert_allclose(fractions, ys, atol=checker.ROUND, rtol=0)
            np.testing.assert_allclose(temperature, 300, atol=1e-10, rtol=0)
            np.testing.assert_allclose(h, 1, atol=checker.ROUND, rtol=0)
            for component in range(state.shape[-1]):
                bad = state.copy()
                bad[0, 3, component] += checker.field_scales(state)[component]*1e-4
                with self.assertRaises(AssertionError):
                    checker.compare("corrupt component", bad, state, state, checker.ROUND)
            bad = state.copy(); bad[0, 0, 0] = -1
            with self.assertRaises(AssertionError):
                checker.decode(name, bad)

    def test_analytic_pulse_and_signed_scalar_transport(self):
        q = checker.pulse_cell_means(160)
        self.assertAlmostEqual(float(q.sum()*40/160), 3.0, places=8)
        x = (np.arange(160)+0.5)*40/160
        for sign in (-1, 1):
            rhs = checker.scalar_rhs(q, 0.25, sign, 3)
            self.assertAlmostEqual(float(rhs.sum()), 0, places=12)
            self.assertAlmostEqual(float(np.sum(x*rhs)/q.sum()), sign, places=10)
        np.testing.assert_allclose(checker.pulse_cell_means(160, 0.1),
                                   checker.pulse_cell_means(160, -0.1)[::-1], atol=1e-13, rtol=0)

    def test_all_step_budget_face_and_admissibility_corruptions(self):
        name = "gas-liquid"
        ys, _, _ = checker.composition(name, np.ones((1, 8))*0.5)
        initial = checker.conservative(name, ys, np.ones((1, 8)))
        nv = initial.shape[-1]
        trace = np.zeros((1, 8, nv+2))
        trace[..., 0] = 1; trace[..., 3] = 300; trace[..., 4:nv] = ys
        stats = np.zeros((2, 15))
        stats[:, :2] = (0.01, 0.1); stats[:, 2:5] = 1
        stats[:, 12] = 300; stats[:, 13:15] = 0.05
        fields = {"statistics": stats, "composition": np.zeros((2, 2*(nv+1)+5)),
                  "q0": initial, "q": initial, "raw_final": initial,
                  "known": np.broadcast_to(initial, (3, *initial.shape)),
                  "solved": np.broadcast_to(initial, (3, *initial.shape))}
        for key in ("W", "E", "S", "N"):
            fields[key] = trace.copy()
        for key in ("L", "R"):
            fields[key] = np.broadcast_to(initial[:, :1], (1, 9, nv)).copy()
        for key in ("B", "T"):
            fields[key] = np.broadcast_to(initial, (2, 8, nv)).copy()
        checker.validate(name, fields, initial, 2)
        with self.assertRaises(AssertionError):
            checker.validate(name, fields, initial, 3)
        for column in (0, 3, 4, 5, 6, 7, nv, nv+1, 2*(nv+1), 2*(nv+1)+1,
                       2*(nv+1)+3, 2*(nv+1)+4):
            bad = {key: value.copy() for key, value in fields.items()}
            bad["composition"][0, column] = 1e5
            with self.assertRaises(AssertionError):
                checker.validate(name, bad, initial, 2)
        for column, value in ((2, -1), (6, 1), (12, 200), (14, -1)):
            bad = {key: value.copy() for key, value in fields.items()}
            bad["statistics"][0, column] = value
            with self.assertRaises(AssertionError):
                checker.validate(name, bad, initial, 2)
        bad = {key: value.copy() for key, value in fields.items()}
        bad["L"][0, 2, 3] *= 1.01
        with self.assertRaises(AssertionError):
            checker.validate(name, bad, initial, 2)


if __name__ == "__main__":
    unittest.main()
