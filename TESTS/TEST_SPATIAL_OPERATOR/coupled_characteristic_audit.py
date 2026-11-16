"""Independent principal-symbol audit of preliminary two-layer face bounds.

This freezes the density/temperature/composition coefficients and linearizes
the existing hydrostatic and interlayer pressure paths on flat terrain.
It samples locally frozen density, layer depths, gravity metric and shear.
These variations are independent local states, NOT gradients of composition or
variable pressure coefficients within a Riemann problem. It does NOT certify
the full variable-composition, wet/dry or curved-bed system, and must not
enable the two-layer production driver.
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


def provisional_bounds(h1, h2, u1, u2, rho1, rho2, rho_air=1.2, g=9.81, metric=1.0):
    """Reproduce the current total-depth envelope, without modifying it."""
    total = h1 + h2
    speeds = [
        (u, np.sqrt(max(0.0, metric * g * (rho - rho_air) / rho * total)))
        for h, u, rho in ((h1, u1, rho1), (h2, u2, rho2)) if h > 0.0
    ]
    if not speeds:
        return 0.0, 0.0
    return min([0.0] + [u - c for u, c in speeds]), max([0.0] + [u + c for u, c in speeds])


def audit_case(label, h1, h2, u1, u2, rho1, rho2, expect_hyperbolic):
    matrix = principal_matrix(h1, h2, u1, u2, rho1, rho2)
    roots = np.linalg.eigvals(matrix)
    scale = max(1.0, np.max(np.abs(roots)))
    imaginary = float(np.max(np.abs(roots.imag)))
    hyperbolic = imaginary <= 1.0e-10 * scale
    if hyperbolic != expect_hyperbolic:
        raise AssertionError(f"{label}: unexpected hyperbolicity ({roots})")

    lower, upper = provisional_bounds(h1, h2, u1, u2, rho1, rho2)
    if hyperbolic:
        tolerance = 1.0e-10 * scale
        if min(roots.real) < lower - tolerance or max(roots.real) > upper + tolerance:
            raise AssertionError(f"{label}: provisional envelope misses real eigenvalues: {roots}")
        if label == "identical_material":
            expected = np.sqrt(9.81 * (rho1 - 1.2) / rho1 * (h1 + h2))
            if not np.allclose(sorted(roots.real), [-expected, 0.0, 0.0, expected], rtol=0, atol=tolerance):
                raise AssertionError("identical layers must recover the barotropic external mode")
    else:
        # CFL bounds on real wave speeds cannot guarantee stability when
        # the principal symbol has non-real eigenvalues.
        if imaginary <= 1.0e-6:
            raise AssertionError(f"{label}: weak or ambiguous complex-mode fixture")
    print(f"COUPLED_CHARACTERISTICS {label}: bounds=({lower:.6g},{upper:.6g}) "
          f"max_imag={imaginary:.6g}", flush=True)



def audit_frozen_parameter_families():
    """Sweep local frozen states and identify nonhyperbolic shear intervals.

    The interval endpoints are grid-resolved diagnostics, not analytic critical
    thresholds. Never classify complex eigenvalues as a CFL-bound violation.
    """
    shears = np.linspace(0.0, 12.0, 241)
    families = (
        ("equal_density", 1000.0, 1000.0),
        ("stable_density", 1800.0, 1400.0),
        ("strong_density_contrast", 1200.0, 600.0),
        ("near_ambient_upper", 1000.0, 2.0),
    )
    total_hyperbolic = 0
    total_nonhyperbolic = 0
    total_bound_failures = 0

    for label, rho1, rho2 in families:
        for h1, h2 in ((1.4, 0.6), (0.3, 1.7)):
            for metric in (0.6, 1.0):
                nonhyperbolic_shears = []
                hyperbolic_count = 0
                misses = 0
                for shear in shears:
                    u1, u2 = 0.5 * shear, -0.5 * shear
                    roots = np.linalg.eigvals(principal_matrix(
                        h1, h2, u1, u2, rho1, rho2, metric=metric))
                    scale = max(1.0, float(np.max(np.abs(roots))))
                    if np.max(np.abs(roots.imag)) > 1.0e-9 * scale:
                        nonhyperbolic_shears.append(shear)
                        continue
                    hyperbolic_count += 1
                    lower, upper = provisional_bounds(
                        h1, h2, u1, u2, rho1, rho2, metric=metric)
                    if (np.min(roots.real) < lower - 1.0e-9 * scale or
                            np.max(roots.real) > upper + 1.0e-9 * scale):
                        misses += 1

                if misses:
                    raise AssertionError(
                        f"{label}: {misses} hyperbolic states outside the wave envelope")
                if not hyperbolic_count:
                    raise AssertionError(f"{label}: no hyperbolic state sampled")
                if not nonhyperbolic_shears:
                    raise AssertionError(f"{label}: missing nonhyperbolic shear window")
                if label == "equal_density" and nonhyperbolic_shears[0] > shears[1]:
                    raise AssertionError("Equal-density shear did not lose hyperbolicity near rest")

                total_hyperbolic += hyperbolic_count
                total_nonhyperbolic += len(nonhyperbolic_shears)
                total_bound_failures += misses
                print(f"COUPLED_SHEAR_SCAN {label} h1={h1:g} h2={h2:g} G={metric:g} "
                      f"hyperbolic={hyperbolic_count} complex={len(nonhyperbolic_shears)} "
                      f"first_complex={nonhyperbolic_shears[0]:g} "
                      f"last_complex={nonhyperbolic_shears[-1]:g}", flush=True)

    print(f"COUPLED_SHEAR_TOTAL hyperbolic={total_hyperbolic} "
          f"complex={total_nonhyperbolic} bound_misses={total_bound_failures}", flush=True)
    if not total_hyperbolic or not total_nonhyperbolic or total_bound_failures:
        raise AssertionError("Frozen-parameter characteristic audit is incomplete")


def run():
    audit_case("identical_material", 1.4, 0.6, 0.0, 0.0, 1000.0, 1000.0, True)
    audit_case("stable_stratification", 1.4, 0.6, 0.0, 0.0, 1800.0, 1400.0, True)
    audit_case("moderate_shear", 1.4, 0.6, 0.5, -0.5, 1800.0, 1400.0, True)
    audit_case("nonhyperbolic_shear", 1.4, 0.6, 1.5, -1.5, 1800.0, 1400.0, False)
    audit_case("large_hyperbolic_shear", 1.4, 0.6, 5.0, -5.0, 1800.0, 1400.0, True)
    audit_frozen_parameter_families()
    print("PASS: principal-symbol bounds, frozen-parameter sweep and shear warning", flush=True)


if __name__ == "__main__":
    run()
