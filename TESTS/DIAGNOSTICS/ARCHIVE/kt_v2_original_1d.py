#!/usr/bin/env python3
"""1-D Kurganov--Petrova (2007) central-upwind shallow-water test.

This is a compact diagnostic implementation of the paper's well-balanced,
positivity-preserving 1-D scheme, applied to the inclined-excavation profile.
The bed is represented by its values at cell vertices and a continuous
piecewise-linear interpolant; the evolved finite-volume variables are
``(w, hu)``, where ``w = h + B``.  The special dry-state reconstruction,
face flux, source quadrature, and SSP-RK3 time integrator follow the paper.

This remains a reduced single-phase model, not a complete solver port.
"""

from dataclasses import dataclass
import argparse
import numpy as np


@dataclass(frozen=True)
class Config:
    nx: int = 120
    dx: float = 0.25
    x0: float = 0.0
    slope_deg: float = 10.0
    excavation_left: float = 10.0
    excavation_right: float = 20.0
    excavation_depth: float = 2.0
    dam_break_x: float = 10.0
    dam_break_depth: float = 1.0
    gravity: float = 9.81
    theta: float = 1.3
    cfl: float = 0.14
    dt0: float = 1.0e-4
    max_dt: float = 2.0e-2
    dt_output: float = 1.0e-1
    dry_tolerance: float = 1.0e-12
    reconstruction: str = 'kp'
    momentum_reconstruction: str = 'hu'


def minmod(*values):
    arrays = np.broadcast_arrays(*[np.asarray(v, dtype=float) for v in values])
    same_sign = np.logical_and.reduce([a > 0.0 for a in arrays]) | np.logical_and.reduce(
        [a < 0.0 for a in arrays]
    )
    magnitudes = np.minimum.reduce([np.abs(a) for a in arrays])
    sign = np.sign(arrays[0])
    return np.where(same_sign, sign * magnitudes, 0.0)


def limited_slope(values, dx, theta):
    """Generalized-minmod reconstruction used in Kurganov & Petrova (2.12)."""
    v = np.asarray(values, dtype=float)
    slope = np.zeros_like(v)
    if v.size >= 3:
        backward = (v[1:-1] - v[:-2]) / dx
        forward = (v[2:] - v[1:-1]) / dx
        centered = (v[2:] - v[:-2]) / (2.0 * dx)
        slope[1:-1] = minmod(theta * backward, centered, theta * forward)
        # Transmissive boundary closure: extend the adjacent one-sided trend.
        slope[0] = (v[1] - v[0]) / dx
        slope[-1] = (v[-1] - v[-2]) / dx
    elif v.size == 2:
        slope[:] = (v[1] - v[0]) / dx
    return slope


def positivity_limited_eta_slope(slope, w, bed_vertices, half_dx):
    """Apply the K&P face-bed positivity correction to one eta slope field."""
    result = np.asarray(slope, dtype=float).copy()
    east_w = w + half_dx * result
    need_east = east_w < bed_vertices[1:]
    result[need_east] = (
        bed_vertices[1:][need_east] - w[need_east]
    ) / half_dx
    west_w = w - half_dx * result
    need_west = west_w < bed_vertices[:-1]
    result[need_west] = (
        w[need_west] - bed_vertices[:-1][need_west]
    ) / half_dx
    return result


def make_bed(cfg):
    """Sample the analytic bed at vertices; use its linear cell averages."""
    xv = cfg.x0 + np.arange(cfg.nx + 1) * cfg.dx
    tan_slope = np.tan(np.deg2rad(cfg.slope_deg))
    weight = np.clip(np.minimum(
        (xv - cfg.excavation_left) / cfg.dx,
        (cfg.excavation_right - xv) / cfg.dx,
    ), 0.0, 1.0)
    bed_vertices = 10.0 - tan_slope * xv - cfg.excavation_depth * weight
    bed_centers = 0.5 * (bed_vertices[:-1] + bed_vertices[1:])
    x = 0.5 * (xv[:-1] + xv[1:])
    bed_slope = (bed_vertices[1:] - bed_vertices[:-1]) / cfg.dx
    return x, xv, bed_vertices, bed_centers, bed_slope


def make_state(cfg, case):
    x, xv, bed_vertices, bed_centers, bed_slope = make_bed(cfg)
    if case == 'excavation':
        # Cell averages of the piecewise-linear initial thickness profile.
        weight_v = np.clip(np.minimum(
            (xv - cfg.excavation_left) / cfg.dx,
            (cfg.excavation_right - xv) / cfg.dx,
        ), 0.0, 1.0)
        h = cfg.excavation_depth * 0.5 * (weight_v[:-1] + weight_v[1:])
    elif case == 'lake':
        # Exact cell-average h for a constant free surface over the linear bed.
        eta0 = np.max(bed_vertices) + 1.0
        h = np.maximum(0.0, eta0 - bed_centers)
    elif case == 'ritter':
        # Ritter's dam-break into a dry bed: flat bottom and a discontinuity
        # aligned with a cell face (the overlap formula also handles offsets).
        bed_vertices = np.zeros_like(bed_vertices)
        bed_centers = np.zeros_like(bed_centers)
        bed_slope = np.zeros_like(bed_slope)
        left_edge = xv[:-1]
        wet_fraction = np.clip((cfg.dam_break_x-left_edge)/cfg.dx, 0.0, 1.0)
        h = cfg.dam_break_depth * wet_fraction
    else:
        raise ValueError(case)

    # K&P evolve cell averages of U=(w, hu), with w=h+B.
    U = np.zeros((cfg.nx, 2))
    U[:, 0] = h + bed_centers
    return x, xv, bed_vertices, bed_centers, bed_slope, U


def ritter_exact_cell_averages(xv, t, cfg):
    """Cell averages of the dry-bed Ritter solution by Gaussian quadrature."""
    xg, wg = np.polynomial.legendre.leggauss(12)
    left = xv[:-1]
    right = xv[1:]
    points = 0.5*(right-left)[:, None]*xg[None, :] + 0.5*(right+left)[:, None]
    c0 = np.sqrt(cfg.gravity*cfg.dam_break_depth)
    if t <= 0.0:
        return cfg.dam_break_depth*np.clip(
            (cfg.dam_break_x-left)/cfg.dx, 0.0, 1.0)
    xi = (points-cfg.dam_break_x)/t
    h = np.zeros_like(points)
    h[points <= cfg.dam_break_x-c0*t] = cfg.dam_break_depth
    fan = ((points > cfg.dam_break_x-c0*t)
           & (points < cfg.dam_break_x+2.0*c0*t))
    c = (2.0*c0-xi[fan])/3.0
    h[fan] = c*c/cfg.gravity
    return 0.5*np.sum(h*wg[None, :], axis=1)


def reconstruct_faces(U, bed_vertices, cfg):
    """Reconstruct (w,hu), applying K&P's conservative dry-cell correction."""
    slopes = np.column_stack((
        limited_slope(U[:, 0], cfg.dx, cfg.theta),
        limited_slope(U[:, 1], cfg.dx, cfg.theta),
    ))
    slope_pick = (0, 0, 0, 0)
    h_cell = U[:, 0] - 0.5 * (bed_vertices[:-1] + bed_vertices[1:])
    h_slope = limited_slope(h_cell, cfg.dx, cfg.theta)
    bed_slope = (bed_vertices[1:] - bed_vertices[:-1]) / cfg.dx
    eta_slope = slopes[:, 0].copy()
    eta2_slope = h_slope + bed_slope
    if cfg.reconstruction in ('minmod_h_eta', 'h_x'):
        # Express the direct h reconstruction as a free-surface slope via
        # eta2_x = h_x + B_x, using the piecewise-linear vertex bed.
        if cfg.reconstruction == 'minmod_h_eta':
            same_sign = ((eta_slope > 0.0) & (eta2_slope > 0.0)) | (
                (eta_slope < 0.0) & (eta2_slope < 0.0)
            )
            eta_pick = same_sign & (np.abs(eta_slope) < np.abs(eta2_slope))
            eta2_pick = same_sign & (np.abs(eta2_slope) < np.abs(eta_slope))
            tie = same_sign & (np.abs(eta_slope) == np.abs(eta2_slope))
            zero = ~(eta_pick | eta2_pick | tie)
            slope_pick = (int(np.count_nonzero(eta_pick)),
                          int(np.count_nonzero(eta2_pick)),
                          int(np.count_nonzero(zero)),
                          int(np.count_nonzero(tie)))
            slopes[:, 0] = minmod(eta_slope, eta2_slope)
        else:
            # h_x-only control: eta2_x=h_x+B_x is the reconstructed eta slope.
            slopes[:, 0] = eta2_slope
    elif cfg.reconstruction not in ('kp', 'face_sign_blend',
                                    'cellwise_h_eta'):
        raise ValueError(f"unknown reconstruction: {cfg.reconstruction}")
    half_dx = 0.5 * cfg.dx
    w = U[:, 0]

    face_eta_weight = np.zeros(len(bed_vertices), dtype=float)
    if cfg.reconstruction == 'cellwise_h_eta':
        # Blend one slope per cell: alpha=1 is eta2_x=h_x+B_x, alpha=0 is
        # eta_x. First make both endpoint reconstructions admissible, then
        # constrain alpha using only the cell's two face traces. The trigger
        # is a candidate h_x diffusive jump directed uphill (and absent from
        # the eta reconstruction); when that jump would wet an eta-dry trace,
        # cap the h_x weight by the available eta-side depth. The minimum of
        # both face caps is used for the whole cell, preserving one linear
        # reconstruction on each cell.
        eta_slope = positivity_limited_eta_slope(
            eta_slope, w, bed_vertices, half_dx
        )
        eta2_slope = positivity_limited_eta_slope(
            eta2_slope, w, bed_vertices, half_dx
        )
        eta_west_h = w - half_dx * eta_slope - bed_vertices[:-1]
        eta_east_h = w + half_dx * eta_slope - bed_vertices[1:]
        h_west_h = w - half_dx * eta2_slope - bed_vertices[:-1]
        h_east_h = w + half_dx * eta2_slope - bed_vertices[1:]
        eta_west_h = np.maximum(eta_west_h, 0.0)
        eta_east_h = np.maximum(eta_east_h, 0.0)
        h_west_h = np.maximum(h_west_h, 0.0)
        h_east_h = np.maximum(h_east_h, 0.0)

        h_weight = np.ones(cfg.nx, dtype=float)
        if cfg.nx > 1:
            bed_centers_local = 0.5 * (bed_vertices[:-1] + bed_vertices[1:])
            bed_face_slope = (bed_centers_local[1:] - bed_centers_local[:-1]) / cfg.dx
            eta_jump = eta_west_h[1:] - eta_east_h[:-1]
            hx_jump = h_west_h[1:] - h_east_h[:-1]
            scale = np.maximum(1.0, np.maximum(
                np.maximum(eta_west_h[1:], eta_east_h[:-1]),
                np.maximum(h_west_h[1:], h_east_h[:-1])))
            tol = 64.0 * np.finfo(float).eps * scale
            hx_diffuses_uphill = hx_jump * bed_face_slope < -tol
            eta_does_not_diffuse_uphill = eta_jump * bed_face_slope >= -tol
            protect = hx_diffuses_uphill & eta_does_not_diffuse_uphill

            # Each cell gets one bound from each of its faces. If eta's trace
            # is dry, the strict non-wetting bound is alpha=0. For a positive
            # eta trace, allow a gradual transition, capped so the blended
            # trace cannot exceed twice the eta-based trace.
            for face in np.flatnonzero(protect):
                left_cell = face
                right_cell = face + 1
                for cell, eta_h, hx_h in (
                        (left_cell, eta_east_h[left_cell], h_east_h[left_cell]),
                        (right_cell, eta_west_h[right_cell], h_west_h[right_cell])):
                    if hx_h <= eta_h + tol[face]:
                        continue
                    if eta_h <= cfg.dry_tolerance:
                        cap = 0.0
                    else:
                        cap = np.clip(eta_h / hx_h, 0.0, 1.0)
                    h_weight[cell] = min(h_weight[cell], float(cap))

        # alpha is the h_x weight; keep it identical at both faces of a cell.
        slopes[:, 0] = eta_slope + h_weight * (eta2_slope - eta_slope)
        face_eta_weight[:-1] = 1.0 - h_weight
        face_eta_weight[-1] = face_eta_weight[-2] if cfg.nx else 1.0
        n_eta = int(np.count_nonzero(h_weight <= 1.0e-14))
        n_blended = int(np.count_nonzero(
            (h_weight > 1.0e-14) & (h_weight < 1.0 - 1.0e-14)))
        n_h = int(np.count_nonzero(h_weight >= 1.0 - 1.0e-14))
        slope_pick = (n_h, n_blended, n_eta, 0)

    # K&P (2.15)-(2.16): adjust only the free-surface slope when needed so
    # each reconstructed face state stays above the unique vertex bed value.
    slopes[:, 0] = positivity_limited_eta_slope(
        slopes[:, 0], w, bed_vertices, half_dx
    )

    west = U - half_dx * slopes
    east = U + half_dx * slopes

    # Two traces at each face; both use the same bed value B_{j+1/2}.
    w_left = np.r_[west[0, 0], east[:, 0]]
    w_right = np.r_[west[:, 0], east[-1, 0]]
    m_left = np.r_[west[0, 1], east[:, 1]]
    m_right = np.r_[west[:, 1], east[-1, 1]]
    h_left = w_left - bed_vertices
    h_right = w_right - bed_vertices

    if cfg.reconstruction == 'face_sign_blend':
        # K&P/eta candidate, after its bed-positivity correction.
        h_eta_left = np.maximum(h_left, 0.0)
        h_eta_right = np.maximum(h_right, 0.0)

        # Direct h_x traces, evaluated at each cell edge and clipped to the
        # nonnegative admissible set. Both candidates use the same vertex bed
        # at a face when converted back to eta.
        h_west = h_cell - half_dx * h_slope
        h_east = h_cell + half_dx * h_slope
        h_h_left = np.r_[h_west[0], h_east]
        h_h_right = np.r_[h_west, h_east[-1]]
        h_h_left = np.maximum(h_h_left, 0.0)
        h_h_right = np.maximum(h_h_right, 0.0)

        # Choose the smallest eta-weight theta in h(theta)=(1-theta)h_x+
        # theta*h_eta for which the dissipative jump cannot point against
        # the eta jump. The same theta is used on both sides of this face and
        # the resulting h traces feed the entire conservative face state.
        delta_eta_depth = h_eta_right - h_eta_left
        delta_h_direct = h_h_right - h_h_left
        delta_scale = np.maximum(1.0, np.maximum(
            np.maximum(np.abs(h_eta_left), np.abs(h_eta_right)),
            np.maximum(np.abs(h_h_left), np.abs(h_h_right))))
        tol = 64.0 * np.finfo(float).eps * delta_scale
        wrong_way = (delta_eta_depth * delta_h_direct < 0.0)
        flat_eta = np.abs(delta_eta_depth) <= tol
        blend = wrong_way | flat_eta
        denom = delta_eta_depth - delta_h_direct
        valid = blend & (np.abs(denom) > tol)
        face_eta_weight[valid] = np.clip(
            -delta_h_direct[valid] / denom[valid], 0.0, 1.0)
        face_eta_weight[flat_eta] = 1.0

        # A trace that eta_x identifies as dry must not be made wet by h_x:
        # that would permit transport across a rim without free-surface
        # overtopping. This is a local lower bound on the same eta weight.
        dry_eta_left = h_eta_left <= cfg.dry_tolerance
        dry_eta_right = h_eta_right <= cfg.dry_tolerance
        h_x_wets_eta_dry_trace = (
            (dry_eta_left & (h_h_left > cfg.dry_tolerance))
            | (dry_eta_right & (h_h_right > cfg.dry_tolerance)))
        face_eta_weight[h_x_wets_eta_dry_trace] = 1.0

        h_left = ((1.0-face_eta_weight)*h_h_left
                  + face_eta_weight*h_eta_left)
        h_right = ((1.0-face_eta_weight)*h_h_right
                   + face_eta_weight*h_eta_right)
        n_h_only = int(np.count_nonzero(face_eta_weight <= 1.0e-14))
        n_blended = int(np.count_nonzero(
            (face_eta_weight > 1.0e-14) & (face_eta_weight < 1.0-1.0e-14)))
        n_eta_only = int(np.count_nonzero(face_eta_weight >= 1.0-1.0e-14))
        slope_pick = (n_h_only, n_blended, n_eta_only, 0)
    return h_left, h_right, m_left, m_right, slope_pick


def desingularized_velocity(h, hu, cfg):
    """K&P equation (2.17), with epsilon=(dx)^4, and its sign preserved."""
    eps = cfg.dx ** 4
    denom = np.sqrt(h**4 + np.maximum(h**4, eps))
    u = np.zeros_like(h)
    np.divide(np.sqrt(2.0) * h * hu, denom, out=u, where=denom > 0.0)
    return u


def kt_rhs(U, bed_vertices, bed_centers, bed_slope, cfg, diagnostics=False):
    h_cell = U[:, 0] - bed_centers
    if np.any(h_cell < -cfg.dry_tolerance):
        raise FloatingPointError(f"negative cell-average depth: {np.min(h_cell):.6e}")

    h_left, h_right, m_left_raw, m_right_raw, slope_pick = reconstruct_faces(
        U, bed_vertices, cfg
    )
    # The correction guarantees nonnegative traces analytically; max only
    # removes round-off at a nominally dry face.
    h_left = np.maximum(h_left, 0.0)
    h_right = np.maximum(h_right, 0.0)
    if cfg.momentum_reconstruction == 'hu':
        # Reference K&P path: reconstruct hu, recover face velocity from the
        # reconstructed depth and momentum, then form the face momentum.
        u_left = desingularized_velocity(h_left, m_left_raw, cfg)
        u_right = desingularized_velocity(h_right, m_right_raw, cfg)
    elif cfg.momentum_reconstruction == 'u':
        # Experiment: reconstruct the desingularized cell-centred velocity
        # directly, then use the same h traces as the mass equation to form
        # the face momentum hu = h*u.
        u_center = desingularized_velocity(h_cell, U[:, 1], cfg)
        u_slope = limited_slope(u_center, cfg.dx, cfg.theta)
        u_west = u_center - 0.5 * cfg.dx * u_slope
        u_east = u_center + 0.5 * cfg.dx * u_slope
        u_left = np.r_[u_west[0], u_east]
        u_right = np.r_[u_west, u_east[-1]]
        u_left = np.where(h_left > cfg.dry_tolerance, u_left, 0.0)
        u_right = np.where(h_right > cfg.dry_tolerance, u_right, 0.0)
    else:
        raise ValueError(f"unknown momentum reconstruction: {cfg.momentum_reconstruction}")
    m_left = h_left * u_left
    m_right = h_right * u_right

    c_left = np.sqrt(cfg.gravity * h_left)
    c_right = np.sqrt(cfg.gravity * h_right)
    a_minus = np.minimum(np.minimum(u_left - c_left, u_right - c_right), 0.0)
    a_plus = np.maximum(np.maximum(u_left + c_left, u_right + c_right), 0.0)
    den = a_plus - a_minus

    f_left = np.column_stack((m_left,
                              m_left * u_left + 0.5 * cfg.gravity * h_left**2))
    f_right = np.column_stack((m_right,
                               m_right * u_right + 0.5 * cfg.gravity * h_right**2))
    U_left = np.column_stack((h_left + bed_vertices, m_left))
    U_right = np.column_stack((h_right + bed_vertices, m_right))
    flux = np.zeros_like(U_left)
    flux_advective = np.zeros_like(U_left)
    flux_dissipative = np.zeros_like(U_left)
    moving = den > 0.0
    flux_advective[moving] = (
        a_plus[moving, None] * f_left[moving]
        - a_minus[moving, None] * f_right[moving]
    ) / den[moving, None]
    flux_dissipative[moving] = (
        (a_plus[moving] * a_minus[moving] / den[moving])[:, None]
        * (U_right[moving] - U_left[moving])
    )
    flux = flux_advective + flux_dissipative

    rhs = -(flux[1:] - flux[:-1]) / cfg.dx
    # K&P source quadrature (2.10), using the same vertex-defined bed.
    rhs[:, 1] -= cfg.gravity * h_cell * bed_slope
    if diagnostics:
        return rhs, dict(flux=flux, a_minus=a_minus, a_plus=a_plus,
                         flux_advective=flux_advective,
                         flux_dissipative=flux_dissipative,
                         h_left=h_left, h_right=h_right,
                         u_left=u_left, u_right=u_right,
                         slope_pick=slope_pick,
                         source=-cfg.gravity * h_cell * bed_slope)
    return rhs


def choose_dt(U, bed_vertices, bed_centers, bed_slope, cfg, t, t_end,
              dt_old, dt_old_old):
    _, d = kt_rhs(U, bed_vertices, bed_centers, bed_slope, cfg, diagnostics=True)
    speed = max(float(np.max(np.maximum(np.abs(d['a_minus']),
                                       np.abs(d['a_plus'])))), 1.0e-12)
    # K&P positivity bound for Forward Euler/SSP methods is CFL <= 1/2.
    dt = min(cfg.max_dt, cfg.cfl * cfg.dx / speed, t_end - t)
    if cfg.dt_output > 0.0:
        next_output = (np.floor(t / cfg.dt_output + 1.0e-12) + 1.0) * cfg.dt_output
        dt = min(dt, next_output - t)
    return min(dt, 1.1 * 0.5 * (dt_old + dt_old_old))


def run(cfg, case, t_end):
    x, xv, bed_vertices, bed_centers, bed_slope, U = make_state(cfg, case)
    U_initial = U.copy()
    h_initial = U[:, 0] - bed_centers
    track_uphill = case == 'excavation'
    uphill = x < cfg.excavation_left if track_uphill else None
    initial_uphill_mass = (float(np.sum(h_initial[uphill]) * cfg.dx)
                           if track_uphill else None)
    max_uphill_mass = initial_uphill_mass
    cumulative_uphill_flux = 0.0
    cumulative_uphill_advective = 0.0
    cumulative_uphill_dissipative = 0.0
    uphill_face = (int(round((cfg.excavation_left - cfg.x0) / cfg.dx))
                   if track_uphill else None)
    t = 0.0
    dt_old = dt_old_old = cfg.dt0
    nstep = 0
    min_h = float(np.min(h_initial))
    slope_pick_counts = np.zeros(4, dtype=np.int64)

    while t < t_end - 1.0e-14:
        dt = choose_dt(U, bed_vertices, bed_centers, bed_slope, cfg,
                       t, t_end, dt_old, dt_old_old)
        L0, d0 = kt_rhs(U, bed_vertices, bed_centers, bed_slope, cfg,
                        diagnostics=True)
        U1 = U + dt * L0
        L1, d1 = kt_rhs(U1, bed_vertices, bed_centers, bed_slope, cfg,
                        diagnostics=True)
        U2 = 0.75 * U + 0.25 * (U1 + dt * L1)
        L2, d2 = kt_rhs(U2, bed_vertices, bed_centers, bed_slope, cfg,
                        diagnostics=True)
        U_next = (1.0 / 3.0) * U + (2.0 / 3.0) * (U2 + dt * L2)

        h_next = U_next[:, 0] - bed_centers
        min_h = min(min_h, float(np.min(h_next)))
        if np.min(h_next) < -cfg.dry_tolerance:
            i = int(np.argmin(h_next))
            return dict(status='negative thickness', t=t, dt=dt, steps=nstep,
                        cell=i + 1, x=x[i], h=float(h_next[i]), min_h=min_h,
                        mass=float(np.sum(U[:, 0] - bed_centers) * cfg.dx),
                        U=U, xgrid=x, bed=bed_centers,
                        bed_vertices=bed_vertices, bed_slope=bed_slope,
                        last_rhs=L2, last_diag=d2)

        # SSP-RK3 quadrature of the interface flux over the full time step.
        if track_uphill:
            cumulative_uphill_flux += dt * (
                d0['flux'][uphill_face, 0] / 6.0
                + d1['flux'][uphill_face, 0] / 6.0
                + 2.0 * d2['flux'][uphill_face, 0] / 3.0
            )
            cumulative_uphill_advective += dt * (
                d0['flux_advective'][uphill_face, 0] / 6.0
                + d1['flux_advective'][uphill_face, 0] / 6.0
                + 2.0 * d2['flux_advective'][uphill_face, 0] / 3.0
            )
            cumulative_uphill_dissipative += dt * (
                d0['flux_dissipative'][uphill_face, 0] / 6.0
                + d1['flux_dissipative'][uphill_face, 0] / 6.0
                + 2.0 * d2['flux_dissipative'][uphill_face, 0] / 3.0
            )
        U = U_next
        # Count only the three SSP-RK stages for accepted time steps;
        # CFL estimator calls are deliberately excluded.
        if cfg.reconstruction in ('minmod_h_eta', 'face_sign_blend',
                                  'cellwise_h_eta'):
            slope_pick_counts += np.asarray(d0['slope_pick'])
            slope_pick_counts += np.asarray(d1['slope_pick'])
            slope_pick_counts += np.asarray(d2['slope_pick'])
        t += dt
        nstep += 1
        dt_old_old, dt_old = dt_old, dt
        if track_uphill:
            max_uphill_mass = max(
                max_uphill_mass,
                float(np.sum(np.maximum(U[:, 0] - bed_centers, 0.0)[uphill])
                      * cfg.dx),
            )

    h_final = U[:, 0] - bed_centers
    result = dict(status='completed', t=t, dt=dt_old, steps=nstep, min_h=min_h,
                mass=float(np.sum(h_final) * cfg.dx), U=U, xgrid=x,
                bed=bed_centers, bed_vertices=bed_vertices, bed_slope=bed_slope,
                U_initial=U_initial,
                slope_pick_counts=slope_pick_counts)
    if track_uphill:
        result.update(initial_uphill_mass=initial_uphill_mass,
                      max_uphill_mass=max_uphill_mass,
                      cumulative_uphill_flux=cumulative_uphill_flux,
                      cumulative_uphill_advective=cumulative_uphill_advective,
                      cumulative_uphill_dissipative=cumulative_uphill_dissipative)
    if case == 'ritter':
        h_exact = ritter_exact_cell_averages(xv, t, cfg)
        error = h_final - h_exact
        result.update(ritter_l1_error=float(np.sum(np.abs(error)) * cfg.dx),
                      ritter_linf_error=float(np.max(np.abs(error))),
                      ritter_exact_mass=float(np.sum(h_exact) * cfg.dx),
                      ritter_exact_h=h_exact)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('excavation', 'lake', 'ritter'),
                        default='excavation')
    parser.add_argument('--t-end', type=float, default=1.5)
    parser.add_argument('--cfl', type=float, default=0.14)
    parser.add_argument('--reconstruction', choices=('kp', 'minmod_h_eta',
                                                      'h_x', 'face_sign_blend',
                                                      'cellwise_h_eta'),
                        default='kp')
    parser.add_argument('--momentum-reconstruction', choices=('hu', 'u'),
                        default='hu',
                        help="reconstruct hu as in K&P ('hu') or reconstruct cell-centred u and set hu_face=h_face*u_face ('u')")
    args = parser.parse_args()
    cfg = Config(cfl=args.cfl, reconstruction=args.reconstruction,
                 momentum_reconstruction=args.momentum_reconstruction)
    result = run(cfg, args.case, args.t_end)
    print('Kurganov-Petrova (2007) 1-D Python diagnostic')
    print(f"case={args.case}, reconstruction={cfg.reconstruction}, "
          f"momentum_reconstruction={cfg.momentum_reconstruction}, "
          f"status={result['status']}, t={result['t']:.10e} s, "
          f"steps={result['steps']}, min(raw h)={result['min_h']:+.10e} m, "
          f"mass={result['mass']:.12e} m^2")
    if result['status'] == 'completed' and args.case == 'excavation':
        print(f"uphill mass initial/max={result['initial_uphill_mass']:.12e}/"
              f"{result['max_uphill_mass']:.12e} m^2; "
              f"increase={result['max_uphill_mass']-result['initial_uphill_mass']:+.6e} m^2; "
              f"integrated rightward flux at uphill edge="
              f"{result['cumulative_uphill_flux']:+.6e} m^2")
    if result['status'] == 'completed' and args.case == 'ritter':
        print(f"Ritter dry-bed reference: h0={cfg.dam_break_depth:.6g} m, "
              f"dam x={cfg.dam_break_x:.6g} m; "
              f"mass error={result['mass']-result['ritter_exact_mass']:+.6e} m^2, "
              f"L1(cell-average h)={result['ritter_l1_error']:.6e} m^2, "
              f"Linf(cell-average h)={result['ritter_linf_error']:.6e} m")
    if result['status'] == 'completed' and cfg.reconstruction in (
            'minmod_h_eta', 'face_sign_blend', 'cellwise_h_eta'):
        eta_pick, eta2_pick, zero_pick, ties = result['slope_pick_counts']
        n_decisions = int(np.sum(result['slope_pick_counts']))
        if cfg.reconstruction == 'minmod_h_eta':
            print(f"local minmod picks over accepted RK stages/cells: "
                  f"eta_x={eta_pick} ({100.0*eta_pick/n_decisions:.2f}%), "
                  f"eta2_x=h_x+B_x={eta2_pick} "
                  f"({100.0*eta2_pick/n_decisions:.2f}%), zero={zero_pick} "
                  f"({100.0*zero_pick/n_decisions:.2f}%), ties={ties} "
                  f"({100.0*ties/n_decisions:.2f}%); total={n_decisions}")
        elif cfg.reconstruction == 'face_sign_blend':
            print(f"face-sign/non-overtopping blend decisions over accepted RK stages/faces: "
                  f"h_x only={eta_pick} ({100.0*eta_pick/n_decisions:.2f}%), "
                  f"partial eta blend={eta2_pick} "
                  f"({100.0*eta2_pick/n_decisions:.2f}%), "
                  f"eta_x fallback={zero_pick} "
                  f"({100.0*zero_pick/n_decisions:.2f}%); total={n_decisions}")
        else:
            h_only, blended, eta_only, _ = result['slope_pick_counts']
            print(f"cellwise h_x-weight decisions over accepted RK stages/cells: "
                  f"h_x only={h_only} ({100.0*h_only/n_decisions:.2f}%), "
                  f"partial blend={blended} ({100.0*blended/n_decisions:.2f}%), "
                  f"eta_x safety limit={eta_only} "
                  f"({100.0*eta_only/n_decisions:.2f}%); total={n_decisions}")
    if result['status'] == 'completed' and args.case == 'excavation':
        print(f"integrated uphill-edge mass flux (advective,dissipative,total)="
              f"({result['cumulative_uphill_advective']:+.6e}, "
              f"{result['cumulative_uphill_dissipative']:+.6e}, "
              f"{result['cumulative_uphill_flux']:+.6e}) m^2")
    if result['status'] == 'completed' and args.case == 'lake':
        h = result['U'][:, 0] - result['bed']
        u = np.divide(result['U'][:, 1], h, out=np.zeros_like(h),
                      where=h > cfg.dry_tolerance)
        h0 = result['U_initial'][:, 0] - result['bed']
        print(f"lake-at-rest drift: max |delta h|="
              f"{np.max(np.abs(h-h0)):.6e} m, max |u|={np.max(np.abs(u)):.6e} m/s")
    if result['status'] != 'completed':
        print(f"first failure: cell={result['cell']}, x={result['x']:.6f} m, "
              f"h={result['h']:+.12e} m, attempted dt={result['dt']:.8e} s")


if __name__ == '__main__':
    main()
