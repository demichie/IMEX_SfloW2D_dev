#!/usr/bin/env python3
"""Reduced 1-D KT/HP-PCCU pressure-split diagnostic.

The script isolates one-phase shallow-water flow (h, h*u), with no settling,
entrainment, rheology, or implicit terms. It keeps separately reconstructed
h, h*u, and u for the inertial advective flux, but removes hydrostatic pressure
from that flux. Hydrostatic momentum is instead integrated as the HP-PCCU
nonconservative path -g*h*d(eta), using limited eta traces and one shared bed
value at each face; both interface and within-cell paths are included. The
central-upwind dissipative jump uses the eta-based face depths, as in HP-PCCU.

This is an exploratory reduced reference, not a complete Python port of the
multiphase IMEX solver. SSPRK2 is used because the isolated case has no implicit
physics.
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
    gravity: float = 9.81
    theta: float = 1.3
    reconstruction_coefficient: float = 1.0
    cfl: float = 0.14
    dt0: float = 1.0e-4
    max_dt: float = 2.0e-2
    dt_output: float = 1.0e-1
    dry_tolerance: float = 1.0e-12
    bed_reconstruction: str = "shared-face"


def minmod(a, b):
    a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
    return np.where(a * b > 0.0, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)


def v2_generalized_minmod_slope(values, dx, theta):
    """Uniform-grid equivalent of v2 ``limit(..., limiter=3)``."""
    v = np.asarray(values, dtype=float)
    slope = np.zeros_like(v)
    if v.size >= 3:
        backward = (v[1:-1] - v[:-2]) / dx
        forward = (v[2:] - v[1:-1]) / dx
        centered = (v[2:] - v[:-2]) / (2.0 * dx)
        slope[1:-1] = minmod(centered, theta * minmod(forward, backward))
    # The v2 boundary reconstruction is BC-dependent. This diagnostic uses
    # transmissive ghosts, which give zero limited slope at the two endpoints.
    return slope


def reconstruct_cell_values(values, cfg):
    slope = v2_generalized_minmod_slope(
        values, cfg.dx, cfg.theta) * cfg.reconstruction_coefficient
    return values - 0.5 * cfg.dx * slope, values + 0.5 * cfg.dx * slope, slope


def bed_and_initial_state(cfg):
    x = cfg.x0 + (np.arange(cfg.nx) + 0.5) * cfg.dx
    base_slope = np.tan(np.deg2rad(cfg.slope_deg))
    # Cell-centred B, as in v2. The Q1 ramp is one cell wide at either edge.
    w_center = np.clip(np.minimum(
        (x - cfg.excavation_left) / cfg.dx,
        (cfg.excavation_right - x) / cfg.dx), 0.0, 1.0)
    B = 10.0 - base_slope * x - cfg.excavation_depth * w_center

    xf = cfg.x0 + np.arange(cfg.nx + 1) * cfg.dx
    w_face = np.clip(np.minimum(
        (xf - cfg.excavation_left) / cfg.dx,
        (cfg.excavation_right - xf) / cfg.dx), 0.0, 1.0)
    B_face = 10.0 - base_slope * xf - cfg.excavation_depth * w_face
    h = cfg.excavation_depth * 0.5 * (w_face[:-1] + w_face[1:])
    q = np.zeros((cfg.nx, 2), dtype=float)
    q[:, 0] = h
    return x, B, B_face, q


def reconstruct(q, B, cfg):
    h = np.maximum(q[:, 0], 0.0)
    hu = q[:, 1].copy()
    u = np.zeros_like(h)
    np.divide(hu, h, out=u, where=h > cfg.dry_tolerance)
    hL, hR, h_slope = reconstruct_cell_values(h, cfg)
    huL, huR, _ = reconstruct_cell_values(hu, cfg)
    uL, uR, _ = reconstruct_cell_values(u, cfg)
    bL, bR, b_slope = reconstruct_cell_values(B, cfg)

    # Match the v2 reconstruction safeguard for exactly empty interior cells:
    # both reconstructed conservative and physical flow states are set dry.
    dry_cell = h <= 0.0
    hL[dry_cell] = 0.0
    hR[dry_cell] = 0.0
    huL[dry_cell] = 0.0
    huR[dry_cell] = 0.0
    uL[dry_cell] = 0.0
    uR[dry_cell] = 0.0

    hL = np.maximum(hL, 0.0)
    hR = np.maximum(hR, 0.0)
    # Reconstructed q=(h, h*u) and auxiliary u are intentionally separate,
    # just as qp_to_qc/qp_interface and qp(idx_u) are separate in v2.
    qL = np.column_stack((hL, huL))
    qR = np.column_stack((hR, huR))
    return qL, qR, uL, uR, bL, bR, b_slope, hL, hR, h_slope


def pccu_interface_flux(q_adv_left, q_adv_right, u_left, u_right,
                        q_diff_left, q_diff_right, eta_left, eta_right,
                        cfg):
    """HP-PCCU-style inertial central-upwind flux and hydrostatic path term."""
    hL = np.maximum(q_adv_left[:, 0], 0.0)
    hR = np.maximum(q_adv_right[:, 0], 0.0)
    cL = np.sqrt(cfg.gravity * hL)
    cR = np.sqrt(cfg.gravity * hR)
    a_minus = np.minimum(np.minimum(u_left - cL, u_right - cR), 0.0)
    a_plus = np.maximum(np.maximum(u_left + cL, u_right + cR), 0.0)
    den = a_plus - a_minus
    fL = np.zeros_like(q_adv_left)
    fR = np.zeros_like(q_adv_right)
    fL[:, 0] = u_left * hL
    fR[:, 0] = u_right * hR
    fL[:, 1] = u_left * q_adv_left[:, 1]
    fR[:, 1] = u_right * q_adv_right[:, 1]
    flux = np.zeros_like(q_adv_left)
    convective = np.zeros_like(q_adv_left)
    dissipative = np.zeros_like(q_adv_left)
    moving = den > 0.0
    convective[moving] = (
        a_plus[moving, None] * fL[moving] - a_minus[moving, None] * fR[moving]
    ) / den[moving, None]
    dissipative[moving] = (
        (a_plus[moving] * a_minus[moving] / den[moving])[:, None]
        * (q_diff_right[moving] - q_diff_left[moving])
    )
    flux[:] = convective + dissipative

    # Match HP-PCCU's stationary equal-free-surface roundoff safeguard.
    eta_scale = np.maximum(1.0, np.maximum(np.abs(eta_left), np.abs(eta_right)))
    stationary_lake = (
        (u_left == 0.0) & (u_right == 0.0)
        & (np.abs(eta_left - eta_right) <= 8.0 * np.finfo(float).eps * eta_scale)
    )
    dissipative[stationary_lake, 0] = 0.0
    flux[:] = convective + dissipative

    # Do not let the numerical face flux drain a dry donor cell. The sign of
    # the net mass flux determines the donor; if its reconstructed advective
    # depth is dry, suppress both transported mass and momentum at this face.
    # The hydrostatic path contribution is assembled separately below and is
    # intentionally left unchanged.
    dry_donor = (
        ((flux[:, 0] > 0.0) & (hL <= cfg.dry_tolerance))
        | ((flux[:, 0] < 0.0) & (hR <= cfg.dry_tolerance))
    )
    convective[dry_donor, :] = 0.0
    dissipative[dry_donor, :] = 0.0
    flux[dry_donor, :] = 0.0

    # For one phase, constant density and G=1, HP-PCCU integrates
    # -g*h*d(eta) along the straight path between eta-based face states.
    path = np.zeros_like(flux)
    path[:, 1] = -cfg.gravity * 0.5 * (q_diff_left[:, 0] + q_diff_right[:, 0]) \
        * (eta_right - eta_left)
    value_left = flux.copy()
    value_right = flux.copy()
    split = den > 1.0e-14
    # A dry side has no fluid column to receive a hydrostatic path force.
    # Keep the path contribution only in the update of a wet adjacent state.
    left_wet = hL > cfg.dry_tolerance
    right_wet = hR > cfg.dry_tolerance
    left_path = split & left_wet
    right_path = split & right_wet
    value_left[left_path] += (
        a_minus[left_path] / den[left_path])[:, None] * path[left_path]
    value_right[right_path] += (
        a_plus[right_path] / den[right_path])[:, None] * path[right_path]
    return (flux, convective, dissipative, path, value_left, value_right,
            a_minus, a_plus, dry_donor)


def rhs(q, B, B_face, cfg, diagnostics=False):
    qL, qR, uL, uR, bL, bR, b_slope, hL, hR, h_slope = reconstruct(q, B, cfg)
    # Each geometric interface sees the east trace of its left cell and the
    # west trace of its right cell; bed remains side-reconstructed as KT input.
    q_if_L = np.vstack((qL[0], qR))
    q_if_R = np.vstack((qL, qR[-1]))
    u_if_L = np.r_[uL[0], uR]
    u_if_R = np.r_[uL, uR[-1]]
    # HP-PCCU path data: reconstruct eta directly, then subtract one shared B
    # value at each face to obtain its left/right hydrostatic endpoint depths.
    h = np.maximum(q[:, 0], 0.0)
    eta = h + B
    etaL, etaR, _ = reconstruct_cell_values(eta, cfg)
    eta_if_L = np.r_[etaL[0], etaR]
    eta_if_R = np.r_[etaL, etaR[-1]]
    if cfg.bed_reconstruction == "center-limited":
        # Preserve the separate left/right bed traces at each internal face.
        b_if_L = np.r_[bL[0], bR]
        b_if_R = np.r_[bL, bR[-1]]
        h_diff_L = np.maximum(eta_if_L - b_if_L, 0.0)
        h_diff_R = np.maximum(eta_if_R - b_if_R, 0.0)
        h_cell_w = np.maximum(etaL - bL, 0.0)
        h_cell_e = np.maximum(etaR - bR, 0.0)
    else:
        # Keep the original single-valued analytic face bed as a control.
        b_if_L = B_face
        b_if_R = B_face
        h_diff_L = np.maximum(eta_if_L - B_face, 0.0)
        h_diff_R = np.maximum(eta_if_R - B_face, 0.0)
        h_cell_w = np.maximum(etaL - B_face[:-1], 0.0)
        h_cell_e = np.maximum(etaR - B_face[1:], 0.0)

    # Scale the advective conservative traces to eta-based depths for the
    # central-upwind jump, inheriting the independently reconstructed u.
    q_diff_L = np.zeros_like(q_if_L)
    q_diff_R = np.zeros_like(q_if_R)
    q_diff_L[:, 0] = h_diff_L
    q_diff_R[:, 0] = h_diff_R
    wet_adv_L = q_if_L[:, 0] > cfg.dry_tolerance
    wet_adv_R = q_if_R[:, 0] > cfg.dry_tolerance
    q_diff_L[wet_adv_L, 1] = (
        h_diff_L[wet_adv_L] / q_if_L[wet_adv_L, 0] * q_if_L[wet_adv_L, 1])
    q_diff_R[wet_adv_R, 1] = (
        h_diff_R[wet_adv_R] / q_if_R[wet_adv_R, 0] * q_if_R[wet_adv_R, 1])

    (flux, flux_conv, flux_diff, path, value_left, value_right,
     a_minus, a_plus, dry_donor) = pccu_interface_flux(
        q_if_L, q_if_R, u_if_L, u_if_R, q_diff_L, q_diff_R,
        eta_if_L, eta_if_R, cfg)

    # Integrate the same nonconservative -g*h*d(eta) form along the path
    # between each cell's west/east eta traces, as HP-PCCU does in addition
    # to the interface paths, using the selected west/east bed traces.
    cell_path = np.zeros_like(q)
    cell_path[:, 1] = -cfg.gravity * 0.5 * (h_cell_w + h_cell_e) * (etaR - etaL)
    cell_path[h <= cfg.dry_tolerance, 1] = 0.0

    # Cell i receives the right-oriented value at its left face and loses the
    # left-oriented value at its right face. q_t has the opposite sign of the
    # assembled spatial term; the integrated cell path is added accordingly.
    dqdt = (value_right[:-1] - value_left[1:] + cell_path) / cfg.dx
    if diagnostics:
        return dqdt, dict(
            flux=flux, flux_conv=flux_conv, flux_diff=flux_diff,
            path=path, cell_path=cell_path,
            value_left=value_left, value_right=value_right,
            a_minus=a_minus, a_plus=a_plus,
            q_if_L=q_if_L, q_if_R=q_if_R,
            u_if_L=u_if_L, u_if_R=u_if_R,
            h_face_L=q_if_L[:, 0], h_face_R=q_if_R[:, 0],
            B_face=B_face, B_trace_L=b_if_L, B_trace_R=b_if_R,
            B_slope=b_slope,
            eta_face_L=eta_if_L, eta_face_R=eta_if_R,
            h_diff_L=h_diff_L, h_diff_R=h_diff_R,
            dry_donor=dry_donor,
        )
    return dqdt


def raw_stage_min_h(q):
    return float(np.min(q[:, 0]))


def step_ssprk2(q, B, B_face, cfg, dt_cap, uphill_face):
    r0, d0 = rhs(q, B, B_face, cfg, diagnostics=True)
    speed0 = max(float(np.max(np.maximum(np.abs(d0['a_minus']),
                                         np.abs(d0['a_plus'])))), 1.0e-12)
    dt = min(dt_cap, cfg.cfl * cfg.dx / speed0)
    dt_initial = dt
    for retry in range(24):
        q1_raw = q + dt * r0
        min_h1_raw = raw_stage_min_h(q1_raw)
        if min_h1_raw >= -1.0e-12:
            q1 = q1_raw.copy()
            q1[:, 0] = np.maximum(q1[:, 0], 0.0)
            r1, d1 = rhs(q1, B, B_face, cfg, diagnostics=True)
            speed1 = max(float(np.max(np.maximum(np.abs(d1['a_minus']),
                                                 np.abs(d1['a_plus'])))), 1.0e-12)
            if dt <= cfg.cfl * cfg.dx / speed1 * (1.0 + 1.0e-12):
                q2_raw = 0.5 * q + 0.5 * (q1 + dt * r1)
                min_h2_raw = raw_stage_min_h(q2_raw)
                if min_h2_raw >= -1.0e-12:
                    q2 = q2_raw.copy()
                    q2[:, 0] = np.maximum(q2[:, 0], 0.0)
                    # In a conservative 1-D update, the mass change in the
                    # uphill subdomain is the net flux through its two edges.
                    # Return per-edge accumulated convective/diffusive fluxes
                    # with the SSPRK2 half-weight applied.
                    wt = 0.5 * dt
                    edge_conv = wt * np.array([
                        d0['flux_conv'][0, 0] + d1['flux_conv'][0, 0],
                        d0['flux_conv'][uphill_face, 0] + d1['flux_conv'][uphill_face, 0],
                    ])
                    edge_diff = wt * np.array([
                        d0['flux_diff'][0, 0] + d1['flux_diff'][0, 0],
                        d0['flux_diff'][uphill_face, 0] + d1['flux_diff'][uphill_face, 0],
                    ])
                    return q2, dt, (min_h1_raw, min_h2_raw), retry, edge_conv, edge_diff
        dt *= 0.5
        if dt < 1.0e-6 * dt_initial:
            raise RuntimeError(
                'positivity requires dt below 1e-6 of the CFL proposal; '
                'the diagnostic is stopping rather than taking a negligible step')
    raise RuntimeError('pressure-split SSPRK2 could not find a positive step')


def step_fortran_rk2(q, B, B_face, cfg, t, dt_old, dt_old_old,
                     t_end, uphill_face):
    """One Fortran n_RK=2 step for a problem with no implicit/source terms.

    The Fortran explicit tableau is [[0, 0], [1, 0]] with weights [1, 0].
    Thus, when the implicit contribution is absent, its two-stage IMEX update
    reduces to forward Euler: stage 2/final state is q^n + dt*rhs(q^n).
    Timestep selection follows compute_spatial_timestep and the main-loop
    growth cap; unlike step_ssprk2, this deliberately does not retry or clip a
    negative thickness, so it can expose the Fortran failure.
    """
    _, d0 = rhs(q, B, B_face, cfg, diagnostics=True)
    max_speed = max(float(np.max(np.maximum(np.abs(d0['a_minus']),
                                           np.abs(d0['a_plus'])))), 1.0e-12)
    dt_cfl = min(cfg.max_dt, cfg.cfl * cfg.dx / max_speed)
    dt = min(dt_cfl, t_end - t)
    if cfg.dt_output > 0.0:
        next_output = (np.floor(t / cfg.dt_output + 1.0e-12) + 1.0) * cfg.dt_output
        dt = min(dt, next_output - t)
    # Same startup/history cap as IMEX_SfloW2D.f90. The initial histories are
    # both dt0; callers update them only after a successful accepted step.
    dt = min(dt, 1.1 * 0.5 * (dt_old + dt_old_old))

    r0, diag = rhs(q, B, B_face, cfg, diagnostics=True)
    q_next = q + dt * r0
    min_h = raw_stage_min_h(q_next)
    if min_h < 0.0:
        cell = int(np.argmin(q_next[:, 0]))
        raise RuntimeError(
            f'negative thickness at Fortran RK stage 2: h={q_next[cell, 0]:+.12e} m '
            f'in cell={cell + 1} (x={(cell + 0.5) * cfg.dx + cfg.x0:.6f} m), '
            f'dt={dt:.8e} s')

    # Compute the same conservative uphill-flux diagnostic as the SSPRK2
    # runner, now using the single forward-Euler stage weight dt.
    edge_conv = dt * np.array([diag['flux_conv'][0, 0],
                               diag['flux_conv'][uphill_face, 0]])
    edge_diff = dt * np.array([diag['flux_diff'][0, 0],
                               diag['flux_diff'][uphill_face, 0]])
    return q_next, dt, min_h, edge_conv, edge_diff


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--t-end', type=float, default=1.5)
    parser.add_argument('--cfl', type=float, default=0.14)
    parser.add_argument('--integrator', choices=('fortran2', 'ssprk2'),
                        default='fortran2',
                        help='fortran2 matches the n_RK=2 explicit Fortran tableau')
    parser.add_argument('--bed-reconstruction',
                        choices=('shared-face', 'center-limited'),
                        default='shared-face',
                        help='shared analytic face bed (control) or cell-centred B '
                             'reconstructed with the h limiter (v2-style)')
    args = parser.parse_args()
    cfg = Config(cfl=args.cfl, bed_reconstruction=args.bed_reconstruction)
    x, B, B_face, q = bed_and_initial_state(cfg)
    uphill = x < cfg.excavation_left
    initial_mass = float(np.sum(q[:, 0]) * cfg.dx)
    initial_uphill = float(np.sum(q[uphill, 0]) * cfg.dx)
    max_uphill = initial_uphill
    cumulative_uphill_conv = 0.0
    cumulative_uphill_diff = 0.0
    uphill_face = int(round((cfg.excavation_left - cfg.x0) / cfg.dx))
    t = 0.0
    nstep = 0
    min_stage_h = 0.0
    dt_old = cfg.dt0
    dt_old_old = cfg.dt0
    retry_count = 0
    while t < args.t_end - 1.0e-14:
        if args.integrator == 'fortran2':
            try:
                q, dt, stage_h, edge_conv, edge_diff = step_fortran_rk2(
                    q, B, B_face, cfg, t, dt_old, dt_old_old,
                    args.t_end, uphill_face)
            except RuntimeError as error:
                print(f'Stopped at t={t:.8e} s: {error}')
                print(f'accepted_steps={nstep}, positivity_retries=0, '
                      f'last_dt={dt_old:.8e} s')
                return
            dt_old_old, dt_old = dt_old, dt
        else:
            try:
                q, dt, stage_h, retries, edge_conv, edge_diff = step_ssprk2(
                    q, B, B_face, cfg, args.t_end - t, uphill_face)
            except RuntimeError as error:
                print(f'Stopped at t={t:.8e} s: {error}')
                return
            retry_count += retries
        t += dt
        nstep += 1
        # Positive values add mass to the uphill region: left-edge inflow
        # minus right-edge outflow.
        cumulative_uphill_conv += edge_conv[0] - edge_conv[1]
        cumulative_uphill_diff += edge_diff[0] - edge_diff[1]
        if np.isscalar(stage_h):
            min_stage_h = min(min_stage_h, float(stage_h))
        else:
            min_stage_h = min(min_stage_h, *stage_h)
        max_uphill = max(max_uphill, float(np.sum(q[uphill, 0]) * cfg.dx))
    final_mass = float(np.sum(q[:, 0]) * cfg.dx)
    h = q[:, 0]
    u = np.divide(q[:, 1], h, out=np.zeros_like(h), where=h > cfg.dry_tolerance)
    print('KT/HP-PCCU pressure-split 1-D inclined-excavation diagnostic')
    print(f'integrator={args.integrator}, bed={args.bed_reconstruction}, '
          f't={t:.8e} s, steps={nstep}, '
          f'min(raw RK-stage h)={min_stage_h:.8e} m, '
          f'positivity_retries={retry_count}')
    print(f'mass initial/final={initial_mass:.12e}/{final_mass:.12e} m^2; '
          f'delta={final_mass-initial_mass:+.6e} m^2')
    print(f'uphill mass initial/max={initial_uphill:.12e}/{max_uphill:.12e} m^2; '
          f'increase={max_uphill-initial_uphill:+.6e} m^2')
    print(f'uphill mass-change contributions ({args.integrator} integrated): '
          f'convective={cumulative_uphill_conv:+.8e} m^2, '
          f'dissipative={cumulative_uphill_diff:+.8e} m^2, '
          f'sum={cumulative_uphill_conv+cumulative_uphill_diff:+.8e} m^2')
    print(f'max |cell u|={np.max(np.abs(u)):.8e} m/s at cell={np.argmax(np.abs(u))+1}, '
          f'h={h[np.argmax(np.abs(u))]:.8e} m')


if __name__ == '__main__':
    main()
