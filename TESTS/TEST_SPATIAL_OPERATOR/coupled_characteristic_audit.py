"""Independent principal-symbol audit of preliminary two-layer face bounds.

This freezes the density/temperature/composition coefficients and linearizes
the existing hydrostatic and interlayer pressure paths on flat terrain.
It does NOT certify the full variable-composition, wet/dry or curved-bed
two-layer system and must not enable the two-layer production driver.
"""

import numpy as np


def principal_matrix(h1, h2, u1, u2, rho1, rho2, rho_air=1.2, g=9.81):
    """Return the homogeneous 1-D primitive Jacobian for (h1,u1,h2,u2).

    Hydrostatic pressure gives gprime_i = g*(rho_i-rho_air)/rho_i.
    L1 receives the interlayer load Gamma2*h2, Gamma2=rho2*gprime2;
    L2 receives gprime2*dh1 through eta2=h1+h2+B.
    """
    g1 = g * (rho1 - rho_air) / rho1
    g2 = g * (rho2 - rho_air) / rho2
    return np.array([
        [u1, h1, 0.0, 0.0],
        [g1, u1, (rho2 / rho1) * g2, 0.0],
        [0.0, 0.0, u2, h2],
        [g2, 0.0, g2, u2],
    ], dtype=float)


def provisional_bounds(h1, h2, u1, u2, rho1, rho2, rho_air=1.2, g=9.81):
    """Reproduce the current total-depth envelope, without modifying it."""
    total = h1 + h2
    speeds = [
        (u, np.sqrt(max(0.0, g * (rho - rho_air) / rho * total)))
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


def run():
    audit_case("identical_material", 1.4, 0.6, 0.0, 0.0, 1000.0, 1000.0, True)
    audit_case("stable_stratification", 1.4, 0.6, 0.0, 0.0, 1800.0, 1400.0, True)
    audit_case("moderate_shear", 1.4, 0.6, 0.5, -0.5, 1800.0, 1400.0, True)
    audit_case("nonhyperbolic_shear", 1.4, 0.6, 1.5, -1.5, 1800.0, 1400.0, False)
    audit_case("large_hyperbolic_shear", 1.4, 0.6, 5.0, -5.0, 1800.0, 1400.0, True)
    print("PASS: principal-symbol eigenvalues, total-depth envelope and shear warning", flush=True)


if __name__ == "__main__":
    run()
