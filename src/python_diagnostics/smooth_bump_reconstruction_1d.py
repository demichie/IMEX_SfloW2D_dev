#!/usr/bin/env python3
"""Compare h_x, eta_x and minmod reconstructions on steady SWASHES bumps.

The script builds the fully-wet steady solutions of the classical shallow-
water bump cases (subcritical, supercritical, and smooth transcritical),
then evaluates the Kurganov--Petrova piecewise-linear reconstruction and its
semidiscrete finite-volume residual on those exact cell averages. It does not
advance the solution in time; this isolates spatial reconstruction error from
boundary transients and time-integration error.
"""

from dataclasses import dataclass
import argparse
import numpy as np


@dataclass(frozen=True)
class BumpCase:
    name: str
    discharge: float
    energy: float
    branch: str


G = 9.81
X_CREsT = 10.0
B_MAX = 0.2
DOMAIN = (0.0, 25.0)
THETA = 1.3
MODES = ("eta_x", "h_x", "minmod")


def bed(x):
    x = np.asarray(x, dtype=float)
    return np.where(np.abs(x - X_CREsT) <= 2.0,
                    B_MAX - 0.05 * (x - X_CREsT) ** 2, 0.0)


def cases():
    # SWASHES subcritical and supercritical benchmark data. The transcritical
    # smooth case is set exactly critical at the bump crest (Fr=1, B=0.2).
    q_sub, h_sub_right = 4.42, 2.0
    e_sub = h_sub_right + q_sub**2 / (2.0 * G * h_sub_right**2)
    q_super, h_super_left = 10.0, 1.0
    e_super = h_super_left + q_super**2 / (2.0 * G * h_super_left**2)
    q_trans = 1.53
    h_critical = (q_trans**2 / G) ** (1.0 / 3.0)
    e_trans = B_MAX + 1.5 * h_critical
    return (
        BumpCase("subcritical", q_sub, e_sub, "high"),
        BumpCase("supercritical", q_super, e_super, "low"),
        BumpCase("transcritical", q_trans, e_trans, "transcritical"),
    )


def depth_from_energy(x, case):
    """Solve h + q^2/(2 g h^2) + B = E on the specified branch."""
    x = np.asarray(x, dtype=float)
    target = case.energy - bed(x)
    hcrit = (case.discharge**2 / G) ** (1.0 / 3.0)
    emin = 1.5 * hcrit
    if np.any(target < emin - 2.0e-12):
        raise ValueError(f"{case.name}: energy is below critical at some point")
    target = np.maximum(target, emin)
    branch = np.full(x.shape, case.branch, dtype=object)
    if case.branch == "transcritical":
        branch = np.where(x < X_CREsT, "high", "low")

    out = np.empty_like(x)
    for index in np.ndindex(x.shape):
        e = target[index]
        side = branch[index]
        if e - emin <= 5.0e-14 * max(1.0, e):
            out[index] = hcrit
            continue

        def residual(h):
            return h + case.discharge**2 / (2.0 * G * h**2) - e

        if side == "high":
            lo, hi = hcrit, max(2.0 * e, 2.0 * hcrit)
            while residual(hi) < 0.0:
                hi *= 2.0
            for _ in range(64):
                mid = 0.5 * (lo + hi)
                if residual(mid) > 0.0:
                    hi = mid
                else:
                    lo = mid
            out[index] = 0.5 * (lo + hi)
        else:
            lo = max(hcrit * 1.0e-12, 1.0e-14)
            hi = hcrit
            for _ in range(64):
                mid = 0.5 * (lo + hi)
                if residual(mid) > 0.0:
                    lo = mid
                else:
                    hi = mid
            out[index] = 0.5 * (lo + hi)
    return out


def limited_slope(v, dx):
    v = np.asarray(v, dtype=float)
    slope = np.zeros_like(v)
    backward = (v[1:-1] - v[:-2]) / dx
    forward = (v[2:] - v[1:-1]) / dx
    centered = 0.5 * (backward + forward)
    slope[1:-1] = minmod(THETA * backward, centered, THETA * forward)
    # Ghost cells are exact constant extensions outside the bump, while the
    # steady state is spatially constant there, so zero boundary slope is apt.
    return slope


def minmod(*values):
    arrays = np.broadcast_arrays(*[np.asarray(v, dtype=float) for v in values])
    same_positive = np.logical_and.reduce([a > 0.0 for a in arrays])
    same_negative = np.logical_and.reduce([a < 0.0 for a in arrays])
    magnitude = np.minimum.reduce([np.abs(a) for a in arrays])
    return np.where(same_positive, magnitude,
                    np.where(same_negative, -magnitude, 0.0))


def quadrature_cell_averages(x_faces, case, nquad=12):
    nodes, weights = np.polynomial.legendre.leggauss(nquad)
    left, right = x_faces[:-1], x_faces[1:]
    xq = (0.5 * (right - left)[:, None] * nodes[None, :]
          + 0.5 * (right + left)[:, None])
    hq = depth_from_energy(xq, case)
    hcell = 0.5 * np.sum(hq * weights[None, :], axis=1)
    # Piecewise-linear cell bed average, consistent with vertex geometry.
    Bface = bed(x_faces)
    Bcell = 0.5 * (Bface[:-1] + Bface[1:])
    return hcell, Bcell, Bface


def reconstruct(h, eta, Bface, dx, mode):
    sh = limited_slope(h, dx)
    seta = limited_slope(eta, dx)
    # In the continuous piecewise-linear vertex bed representation, B_x is
    # the exact cellwise face-to-face slope.
    sB = (Bface[1:] - Bface[:-1]) / dx
    seta2 = sh + sB
    if mode == "eta_x":
        slope = seta.copy()
    elif mode == "h_x":
        slope = seta2.copy()
    elif mode == "minmod":
        slope = minmod(seta, seta2)
    else:
        raise ValueError(mode)

    half = 0.5 * dx
    # K&P positivity correction against the single continuous bed value at
    # each face. The input cases are fully wet, but this remains important at
    # the steep shoulders of the bump.
    east = eta + half * slope
    bad_east = east < Bface[1:]
    slope[bad_east] = (Bface[1:][bad_east] - eta[bad_east]) / half
    west = eta - half * slope
    bad_west = west < Bface[:-1]
    slope[bad_west] = (eta[bad_west] - Bface[:-1][bad_west]) / half

    eta_w = eta - half * slope
    eta_e = eta + half * slope
    h_w = np.maximum(eta_w - Bface[:-1], 0.0)
    h_e = np.maximum(eta_e - Bface[1:], 0.0)
    return slope, h_w, h_e


def semidiscrete_residual(h, Bcell, Bface, dx, case, mode):
    eta = h + Bcell
    slope, h_w, h_e = reconstruct(h, eta, Bface, dx, mode)
    # Reconstructed endpoint values of w=eta and conservative momentum hu=q.
    eta_w = eta - 0.5 * dx * slope
    eta_e = eta + 0.5 * dx * slope
    m = np.full_like(h, case.discharge)
    # Two ghost cells at each end provide the full limiter stencil. The state
    # is analytically constant outside the bump.
    nx = h.size
    xfaces = np.linspace(DOMAIN[0] - 2 * dx, DOMAIN[1] + 2 * dx, nx + 5)
    xc = 0.5 * (xfaces[:-1] + xfaces[1:])
    hghost = depth_from_energy(xc, case)
    Bfg = bed(xfaces)
    Bcg = 0.5 * (Bfg[:-1] + Bfg[1:])
    # Replace interior ghost-generated states by the quadrature averages used
    # for the FV data; the analytic ghost cells remain exact constant states.
    hfull = hghost.copy()
    hfull[2:-2] = h
    Bcfull = Bcg
    etafull = hfull + Bcfull
    slope_h = limited_slope(hfull, dx)
    slope_eta = limited_slope(etafull, dx)
    slope_B = (Bfg[1:] - Bfg[:-1]) / dx
    slope_eta2 = slope_h + slope_B
    if mode == "eta_x":
        s = slope_eta.copy()
    elif mode == "h_x":
        s = slope_eta2.copy()
    else:
        s = minmod(slope_eta, slope_eta2)
    half = 0.5 * dx
    east = etafull + half * s
    bad = east < Bfg[1:]
    s[bad] = (Bfg[1:][bad] - etafull[bad]) / half
    west = etafull - half * s
    bad = west < Bfg[:-1]
    s[bad] = (etafull[bad] - Bfg[:-1][bad]) / half

    etaL = etafull - half * s
    etaR = etafull + half * s
    hL = np.maximum(etaL - Bfg[:-1], 0.0)
    hR = np.maximum(etaR - Bfg[1:], 0.0)
    # Momentum is constant and hence has zero limited slope.
    mL = mR = np.full_like(hfull, case.discharge)
    hLface = np.r_[hL[0], hR]
    hRface = np.r_[hL, hR[-1]]
    mLface = np.r_[mL[0], mR]
    mRface = np.r_[mL, mR[-1]]
    uL = np.divide(mLface, hLface, out=np.zeros_like(hLface), where=hLface > 0.0)
    uR = np.divide(mRface, hRface, out=np.zeros_like(hRface), where=hRface > 0.0)
    cL, cR = np.sqrt(G * hLface), np.sqrt(G * hRface)
    amin = np.minimum(np.minimum(uL - cL, uR - cR), 0.0)
    amax = np.maximum(np.maximum(uL + cL, uR + cR), 0.0)
    den = amax - amin
    # The conservative state in the mass equation is h (not eta). The bed is
    # handled by the well-balanced momentum source, not by diffusing eta in
    # the mass flux.
    UL = np.column_stack((hLface, mLface))
    UR = np.column_stack((hRface, mRface))
    FL = np.column_stack((mLface, mLface * uL + 0.5 * G * hLface**2))
    FR = np.column_stack((mRface, mRface * uR + 0.5 * G * hRface**2))
    flux = np.zeros_like(UL)
    moving = den > 0.0
    flux[moving] = (
        amax[moving, None] * FL[moving] - amin[moving, None] * FR[moving]
        + (amax[moving] * amin[moving] / den[moving])[:, None]
        * (UR[moving] - UL[moving])
    ) / den[moving, None]

    # Physical cells occupy [2:nx+2] in the extended array.
    idx = np.arange(2, nx + 2)
    rhs = -(flux[idx + 1] - flux[idx]) / dx
    rhs[:, 1] -= G * h * ((Bface[1:] - Bface[:-1]) / dx)
    return rhs, (h_w, h_e, slope)


def run_case(case, nx):
    dx = (DOMAIN[1] - DOMAIN[0]) / nx
    # Two exact constant-state ghost cells on either side.
    xfaces = np.linspace(DOMAIN[0] - 2 * dx, DOMAIN[1] + 2 * dx, nx + 5)
    xcenter = 0.5 * (xfaces[:-1] + xfaces[1:])
    hfull = depth_from_energy(xcenter, case)
    h = hfull[2:-2]
    # Exact quadrature averages in active cells.
    xphys = np.linspace(DOMAIN[0], DOMAIN[1], nx + 1)
    h, Bcell, Bface = quadrature_cell_averages(xphys, case)
    rhs_metrics = {}
    trace_metrics = {}
    xmid = 0.5 * (xphys[:-1] + xphys[1:])
    h_exact_faces = depth_from_energy(xphys, case)
    for mode in MODES:
        rhs, (h_w, h_e, slope) = semidiscrete_residual(
            h, Bcell, Bface, dx, case, mode
        )
        trace_error = np.r_[h_w - h_exact_faces[:-1], h_e - h_exact_faces[1:]]
        trace_ref = np.r_[h_exact_faces[:-1], h_exact_faces[1:]]
        residual_mass = float(np.max(np.abs(rhs[:, 0])))
        residual_mom = float(np.max(np.abs(rhs[:, 1])))
        rhs_metrics[mode] = (residual_mass, residual_mom)
        trace_metrics[mode] = (
            float(np.mean(np.abs(trace_error))),
            float(np.max(np.abs(trace_error))),
            float(np.sqrt(np.mean(trace_error**2))),
        )
    return rhs_metrics, trace_metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, nargs="+", default=[100, 200, 400])
    args = parser.parse_args()
    print("Fully wet steady bump reconstruction study (K&P vertex-bed FV residual)")
    print("Residuals are max |d(eta)/dt| [m/s] and max |d(hu)/dt| [m^2/s^2].")
    print("Trace errors are mean / max / RMS depth error at cell faces [m].")
    for case in cases():
        print(f"\n{case.name}: discharge={case.discharge:.6g}, energy={case.energy:.8g}")
        print(" nx    mode       mass-resid      mom-resid       trace-L1       trace-Linf      trace-RMS")
        for nx in args.nx:
            residual, trace = run_case(case, nx)
            for mode in MODES:
                rm, rp = residual[mode]
                e1, einf, erms = trace[mode]
                print(f"{nx:4d}  {mode:8s}  {rm:12.4e}  {rp:12.4e}  "
                      f"{e1:12.4e}  {einf:12.4e}  {erms:12.4e}")


if __name__ == "__main__":
    main()
