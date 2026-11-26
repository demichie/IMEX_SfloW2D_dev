"""Audit frozen two-layer gravitational eigenvalues and the total-depth CFL bound.

Freeze density, temperature and composition locally and linearize the implemented
hydrostatic and interlayer pressure paths on flat terrain. Sample positive depths,
density contrasts and interlayer shear, without changing solver coefficients.
The envelope is assessed only for *real* eigenvalues under stable density ordering.
A real spectrum alone does not establish strong hyperbolicity (equal-density
states are defective), and no result here certifies variable-composition,
wet/dry or curved-bed dynamics.
"""

import numpy as np


def principal_matrix(h1, h2, u1, u2, rho1, rho2, rho_air=1.2, g=9.81, metric=1.0):
    """Return the homogeneous 1-D primitive Jacobian for (h1,u1,h2,u2).

    Hydrostatic pressure gives gprime_i = g*(rho_i-rho_air)/rho_i.
    L1 receives the interlayer load Gamma2*h2, Gamma2=rho2*gprime2;
    L2 receives gprime2*dh1 through eta2=h1+h2+B.
    """
    g1 = metric * g * (rho1 - rho_air) / rho1
    g2 = metric * g * (rho2 - rho_air) / rho2
    return np.array([
        [u1, h1, 0.0, 0.0],
        [g1, u1, (rho2 / rho1) * g2, 0.0],
        [0.0, 0.0, u2, h2],
        [g2, 0.0, g2, u2],
    ], dtype=float)


def total_depth_bounds(h1, h2, u1, u2, rho1, rho2, rho_air=1.2, g=9.81, metric=1.0):
    """Reproduce the current total-depth envelope, without modifying it."""
    total = h1 + h2
    speeds = [
        (u, np.sqrt(max(0.0, metric * g * (rho - rho_air) / rho * total)))
        for h, u, rho in ((h1, u1, rho1), (h2, u2, rho2)) if h > 0.0
    ]
    if not speeds:
        return 0.0, 0.0
    return min([0.0] + [u - c for u, c in speeds]), max([0.0] + [u + c for u, c in speeds])


def audit_case(label, h1, h2, u1, u2, rho1, rho2, expect_real_spectrum):
    matrix = principal_matrix(h1, h2, u1, u2, rho1, rho2)
    roots = np.linalg.eigvals(matrix)
    scale = max(1.0, np.max(np.abs(roots)))
    imaginary = float(np.max(np.abs(roots.imag)))
    # Coincident internal modes can split by O(sqrt(machine epsilon)).
    has_real_spectrum = imaginary <= 1.0e-7 * scale
    if has_real_spectrum != expect_real_spectrum:
        raise AssertionError(f"{label}: unexpected reality of gravitational spectrum ({roots})")

    lower, upper = total_depth_bounds(h1, h2, u1, u2, rho1, rho2)
    if has_real_spectrum:
        tolerance = 1.0e-10 * scale
        if min(roots.real) < lower - tolerance or max(roots.real) > upper + tolerance:
            raise AssertionError(f"{label}: total-depth envelope misses real eigenvalues: {roots}")
        if label == "identical_material":
            expected = np.sqrt(9.81 * (rho1 - 1.2) / rho1 * (h1 + h2))
            if not np.allclose(sorted(roots.real), [-expected, 0.0, 0.0, expected],
                               rtol=0, atol=1.0e-7 * scale):
                raise AssertionError("identical layers must recover the barotropic external mode")
    else:
        # CFL bounds on real wave speeds cannot guarantee stability when
        # the principal symbol has non-real eigenvalues.
        if imaginary <= 1.0e-6:
            raise AssertionError(f"{label}: weak or ambiguous complex-mode fixture")
    print(f"COUPLED_CHARACTERISTICS {label}: bounds=({lower:.6g},{upper:.6g}) "
          f"max_imag={imaginary:.6g}", flush=True)



def audit_equal_density_defect():
    """Check the double zero-slip internal mode has only one eigenvector."""
    matrix = principal_matrix(1.4, 0.6, 0.0, 0.0, 1000.0, 1000.0)
    roots = np.linalg.eigvals(matrix)
    scale = max(1.0, float(np.max(np.abs(roots))))
    zero_roots = np.count_nonzero(np.abs(roots) < 1.0e-7 * scale)
    if zero_roots != 2 or np.linalg.matrix_rank(matrix, tol=1.0e-10 * scale) != 3:
        raise AssertionError("Equal-density zero-slip internal mode must be defective")
    print("COUPLED_CHARACTERISTICS equal_density_defect: "
          "algebraic_multiplicity=2 geometric_multiplicity=1", flush=True)


def audit_frozen_parameter_families():
    """Sweep local frozen states and identify nonhyperbolic shear intervals.

    The interval endpoints are grid-resolved diagnostics, not analytic critical
    thresholds. Never classify complex eigenvalues as a CFL-bound violation.
    Counting real roots is not a test of eigenvector completeness.
    """
    shears = np.linspace(0.0, 12.0, 241)
    families = (
        ("equal_density", 1000.0, 1000.0),
        ("stable_density", 1800.0, 1400.0),
        ("strong_density_contrast", 1200.0, 600.0),
        ("near_ambient_upper", 1000.0, 2.0),
    )
    total_real = 0
    total_complex = 0
    total_bound_failures = 0

    for label, rho1, rho2 in families:
        for h1, h2 in ((1.4, 0.6), (0.3, 1.7)):
            for metric in (0.6, 1.0):
                complex_shears = []
                real_count = 0
                misses = 0
                for shear in shears:
                    u1, u2 = 0.5 * shear, -0.5 * shear
                    roots = np.linalg.eigvals(principal_matrix(
                        h1, h2, u1, u2, rho1, rho2, metric=metric))
                    scale = max(1.0, float(np.max(np.abs(roots))))
                    # Double roots at zero shear are numerically ill-conditioned.
                    if np.max(np.abs(roots.imag)) > 1.0e-7 * scale:
                        complex_shears.append(shear)
                        continue
                    real_count += 1
                    lower, upper = total_depth_bounds(
                        h1, h2, u1, u2, rho1, rho2, metric=metric)
                    if (np.min(roots.real) < lower - 1.0e-9 * scale or
                            np.max(roots.real) > upper + 1.0e-9 * scale):
                        misses += 1

                if misses:
                    raise AssertionError(
                        f"{label}: {misses} real-spectrum states outside the wave envelope")
                if not real_count:
                    raise AssertionError(f"{label}: no real-spectrum state sampled")
                if not complex_shears:
                    raise AssertionError(f"{label}: missing nonhyperbolic shear window")
                if label == "equal_density" and complex_shears[0] > shears[1]:
                    raise AssertionError("Equal-density shear did not lose hyperbolicity near rest")

                total_real += real_count
                total_complex += len(complex_shears)
                total_bound_failures += misses
                print(f"COUPLED_SHEAR_SCAN {label} h1={h1:g} h2={h2:g} G={metric:g} "
                      f"real_spectrum={real_count} complex={len(complex_shears)} "
                      f"first_complex={complex_shears[0]:g} "
                      f"last_complex={complex_shears[-1]:g}", flush=True)

    print(f"COUPLED_SHEAR_TOTAL real_spectrum={total_real} "
          f"complex={total_complex} bound_misses={total_bound_failures}", flush=True)
    if not total_real or not total_complex or total_bound_failures:
        raise AssertionError("Frozen-parameter characteristic audit is incomplete")


def run():
    audit_case("identical_material", 1.4, 0.6, 0.0, 0.0, 1000.0, 1000.0, True)
    audit_equal_density_defect()
    audit_case("inverted_density", 1.4, 0.6, 0.0, 0.0, 1000.0, 1200.0, False)
    audit_case("stable_stratification", 1.4, 0.6, 0.0, 0.0, 1800.0, 1400.0, True)
    audit_case("moderate_shear", 1.4, 0.6, 0.5, -0.5, 1800.0, 1400.0, True)
    audit_case("nonhyperbolic_shear", 1.4, 0.6, 1.5, -1.5, 1800.0, 1400.0, False)
    audit_case("large_hyperbolic_shear", 1.4, 0.6, 5.0, -5.0, 1800.0, 1400.0, True)
    audit_frozen_parameter_families()
    print("PASS: real-root envelope, defective/inverted density and shear checks", flush=True)


if __name__ == "__main__":
    run()
