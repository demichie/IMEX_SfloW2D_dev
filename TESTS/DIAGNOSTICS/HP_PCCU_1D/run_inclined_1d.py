#!/usr/bin/env python3
"""Run a 1-D inclined-excavation diagnostic with the Python HP-PCCU core."""

from pathlib import Path
import argparse
import sys

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parent))

import hp_pccu_1d_core as core  # noqa: E402


NX = 120
DX = 0.25
X0 = 0.0
SLOPE = np.tan(np.deg2rad(10.0))
RHO = 1000.0
EXCAVATION_LEFT = 10.0
EXCAVATION_RIGHT = 20.0
EXCAVATION_DEPTH = 2.0
DEFAULT_T_END = 1.0

PARAMS = core.SolverParams(
    rho_c=RHO, rho_s=RHO, rho_a=0.0, cfl=0.24,
    limiter=3, theta=1.3, reconstr_coeff=1.0,
    slope_correction=False, curvature_term=False,
    dry_h_tol=1.0e-10, zero_transport_at_rest=True,
).checked()


def make_case():
    """Build a sloping bed with a Q1-ramped excavation and initially filled pit."""
    x_face = X0 + np.arange(NX + 1) * DX
    x_center = X0 + (np.arange(NX) + 0.5) * DX

    # The ramp is one cell wide, matching the Q1 geometry used by the 2-D case.
    face_weight = np.clip(
        np.minimum((x_face - EXCAVATION_LEFT) / DX,
                   (EXCAVATION_RIGHT - x_face) / DX), 0.0, 1.0)
    bed = core.continuous_bed_geometry(
        10.0 - SLOPE * x_face - EXCAVATION_DEPTH * face_weight)

    # Cell-average the nodal Q1 weight, as in the existing N7 excavation case.
    h = EXCAVATION_DEPTH * 0.5 * (face_weight[:-1] + face_weight[1:])
    q = core.state_from_h_alpha_hu(h, np.zeros(NX), np.zeros(NX), PARAMS)
    return bed, q, x_center


def run(scale, t_end, local_slope_blend=False, no_leak_rim_gate=False,
        eta_only_transport_reconstruction=False, integrator="ssprk2",
        zero_velocity_on_dry_trace=False,
        zero_advective_flux_from_dry_donor=False):
    params = core.SolverParams(**{
        **PARAMS.__dict__, "hp_blend_scale": scale,
        "hp_blend_diagnostic": True,
        "local_slope_blend": local_slope_blend,
        "no_leak_rim_gate": no_leak_rim_gate,
        "eta_only_transport_reconstruction": eta_only_transport_reconstruction,
        "zero_velocity_on_dry_trace": zero_velocity_on_dry_trace,
        "zero_advective_flux_from_dry_donor": zero_advective_flux_from_dry_donor,
    }).checked()
    bed, q, x_center = make_case()
    initial_reconstruction = core.assemble_face_states(
        q, bed, 0.0, params, core.transmissive_boundary)
    alpha0 = initial_reconstruction["local_blend_alpha"]
    infeasible0 = initial_reconstruction["local_blend_infeasible"]
    uphill = x_center < EXCAVATION_LEFT
    initial_uphill_mass = float(np.sum(q[uphill, 0]) * DX)
    uphill_face = int(round((EXCAVATION_LEFT-X0)/DX))
    reports = [(0.0, initial_uphill_mass, 0.0, 0.0, 0.0, 0.0)]
    cumulative_uphill_advective = 0.0
    cumulative_uphill_dissipative = 0.0
    cumulative_left_advective = 0.0
    cumulative_left_dissipative = 0.0
    maximum_uphill_excess = 0.0
    minimum_accepted_h = float(np.min(core.raw_thickness(q, params)))
    euler_retries = 0
    report_times = sorted(set((min(0.1, t_end), min(0.5, t_end), t_end)))
    targets = iter(report_times)
    target = next(targets, None)
    t = 0.0
    steps = 0

    while target is not None:
        while t < target - 1.0e-12:
            if integrator == "ssprk2":
                q, diag = core.step_ssprk2(
                    q, bed, DX, t, target - t, params,
                    core.transmissive_boundary, return_diagnostics=True)
            elif integrator == "euler":
                # One forward-Euler update per time step. Reduce dt only if
                # this single explicit update would violate h >= 0.
                cache = core.prepare_stage_face_cache(
                    q, bed, t, params, core.transmissive_boundary,
                    q_is_sanitized=True)
                dt_cfl = core.max_dt(
                    q, bed, DX, t, params, core.transmissive_boundary,
                    stage_cache=cache)
                dt = min(target - t, 0.99 * float(dt_cfl))
                retries = 0
                tol_h = 5.0e-13 * max(1.0, float(np.max(cache["data"]["prim_q"][3])))
                while True:
                    rhs, rhs_diag = core.rhs_spatial(
                        q, bed, DX, t, params, core.transmissive_boundary,
                        return_diagnostics=True, stage_cache=cache)
                    h_start = core.raw_thickness(q, params)
                    h_rate = core.raw_thickness(rhs, params)
                    dry_outflow = (h_start <= 1.0e-14) & (h_rate < 0.0)
                    if np.any(dry_outflow):
                        i_bad = int(np.flatnonzero(dry_outflow)[0])
                        f_left = i_bad
                        f_right = i_bad + 1
                        fluxes = []
                        for face, side in ((f_left, "left"),
                                           (f_right, "right")):
                            data = cache["data"]
                            fluxes.append(
                                f"{side} face x={X0 + face*DX:.6f}: "
                                f"(hL,uL)=({core.raw_thickness(data['q_advL_face'][face:face+1], params)[0]:.3e},"
                                f"{data['u_advL_face'][face]:+.3e}), "
                                f"(hR,uR)=({core.raw_thickness(data['q_advR_face'][face:face+1], params)[0]:.3e},"
                                f"{data['u_advR_face'][face]:+.3e}), "
                                f"Hadv={rhs_diag['H_advective'][face, 0]:+.6e}, "
                                f"Hdiff={rhs_diag['H_dissipative'][face, 0]:+.6e}, "
                                f"H={rhs_diag['H'][face, 0]:+.6e}")
                        raise RuntimeError(
                            "Forward-Euler positivity failure is already present "
                            "in the semidiscrete RHS: "
                            f"t={t:.12g}, cell={i_bad + 1}, "
                            f"x={X0 + (i_bad + 0.5)*DX:.6f} m, "
                            f"h={h_start[i_bad]:.12e} m, "
                            f"dh/dt={h_rate[i_bad]:+.12e} m/s; "
                            + "; ".join(fluxes))
                    q_raw = q + dt * rhs
                    h_raw = core.raw_thickness(q_raw, params)
                    min_h = float(np.min(h_raw))
                    if min_h >= -tol_h:
                        q = core.sanitize(q_raw, params)
                        minimum_accepted_h = min(
                            minimum_accepted_h,
                            float(np.min(core.raw_thickness(q, params))))
                        diag = {
                            "dt": dt,
                            "integrated_flux_advective": dt * rhs_diag["H_advective"][:, 0],
                            "integrated_flux_dissipative": dt * rhs_diag["H_dissipative"][:, 0],
                        }
                        euler_retries += retries
                        break
                    retries += 1
                    if retries > 30 or dt <= 1.0e-15:
                        i_min = int(np.argmin(h_raw))
                        raise RuntimeError(
                            "Forward-Euler positivity check failed: "
                            f"t={t:.12g}, dt={dt:.12g}, min_h={min_h:.12e} "
                            f"at cell={i_min + 1}, "
                            f"x={X0 + (i_min + 0.5)*DX:.6f} m, "
                            f"after {retries} retries")
                    dt *= 0.5
            # The left/uphill subdomain gains mass from the x=10 m interface
            # when its +x-oriented numerical flux is negative.
            cumulative_uphill_advective -= float(
                diag["integrated_flux_advective"][uphill_face])
            cumulative_uphill_dissipative -= float(
                diag["integrated_flux_dissipative"][uphill_face])
            cumulative_left_advective += float(
                diag["integrated_flux_advective"][0])
            cumulative_left_dissipative += float(
                diag["integrated_flux_dissipative"][0])
            t += diag["dt"]
            steps += 1
        uphill_mass = float(np.sum(q[uphill, 0]) * DX)
        maximum_uphill_excess = max(
            maximum_uphill_excess, uphill_mass - initial_uphill_mass)
        reports.append((target, uphill_mass, cumulative_uphill_advective,
                        cumulative_uphill_dissipative, cumulative_left_advective,
                        cumulative_left_dissipative))
        target = next(targets, None)

    total_mass = float(np.sum(q[:, 0]) * DX)
    return (reports, maximum_uphill_excess, steps, total_mass, alpha0,
            infeasible0, minimum_accepted_h, euler_retries)


def trial_raw_euler(scale, dt, zero_velocity_on_dry_trace=False,
                    zero_advective_flux_from_dry_donor=False):
    """Inspect one un-clipped Euler stage before the safe SSPRK2 wrapper."""
    params = core.SolverParams(**{
        **PARAMS.__dict__, "hp_blend_scale": scale,
        "hp_blend_diagnostic": True,
        "zero_velocity_on_dry_trace": zero_velocity_on_dry_trace,
        "zero_advective_flux_from_dry_donor": zero_advective_flux_from_dry_donor,
    }).checked()
    bed, q, _ = make_case()
    cache = core.prepare_stage_face_cache(q, bed, 0.0, params,
                                          core.transmissive_boundary)
    rhs = core.rhs_spatial(q, bed, DX, 0.0, params,
                           core.transmissive_boundary, stage_cache=cache)
    trial = q + dt * rhs
    h_raw = core.raw_thickness(trial, params)
    i_min = int(np.argmin(h_raw))
    print(f"  raw Euler trial dt={dt:.8e} s: min_h={h_raw[i_min]:+.12e} m "
          f"at cell={i_min + 1}, x={X0 + (i_min + 0.5)*DX:.6f} m")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t-end", type=float, default=DEFAULT_T_END)
    parser.add_argument("--scale", type=float, default=1.0,
                        help="diagnostic HP blend scale (default: 1)")
    parser.add_argument("--reconstruction",
                        choices=("baseline", "local-slope-blend", "eta-only"),
                        default="baseline",
                        help="select baseline, experimental h/eta slope blend, or eta_x-only traces")
    parser.add_argument("--integrator", choices=("ssprk2", "euler"),
                        default="ssprk2",
                        help="time integrator; euler performs one explicit update per step")
    parser.add_argument("--zero-velocity-on-dry-trace", action="store_true",
                        help="set each advective face velocity to zero when its own h trace is dry")
    parser.add_argument("--zero-advective-flux-from-dry-donor", action="store_true",
                        help="zero the full advective flux vector if its mass-flux donor trace is dry")
    parser.add_argument("--no-leak-rim-gate", action="store_true",
                        help="apply a face-local dry-rim/no-overtopping gate")
    parser.add_argument("--trial-dt", type=float, default=8.6369206229238125e-3,
                        help="report the raw first-stage Euler thickness at this dt")
    args = parser.parse_args()
    if args.t_end <= 0.0:
        parser.error("--t-end must be positive")
    if not 0.0 <= args.scale <= 1.0:
        parser.error("--scale must be between zero and one")

    local_blend = args.reconstruction == "local-slope-blend"
    eta_only = args.reconstruction == "eta-only"
    (reports, max_uphill, steps, total_mass, alpha0, infeasible0,
     min_h, euler_retries) = run(
        args.scale, args.t_end, local_slope_blend=local_blend,
        no_leak_rim_gate=args.no_leak_rim_gate,
        eta_only_transport_reconstruction=eta_only,
        integrator=args.integrator,
        zero_velocity_on_dry_trace=args.zero_velocity_on_dry_trace,
        zero_advective_flux_from_dry_donor=args.zero_advective_flux_from_dry_donor)
    trial_raw_euler(args.scale, args.trial_dt,
                    zero_velocity_on_dry_trace=args.zero_velocity_on_dry_trace,
                    zero_advective_flux_from_dry_donor=args.zero_advective_flux_from_dry_donor)
    initial_uphill_mass = reports[0][1]
    print(f"1-D, integrator={args.integrator}, reconstruction={args.reconstruction}, "
          f"zero_velocity_on_dry_trace={args.zero_velocity_on_dry_trace}, "
          f"zero_advective_flux_from_dry_donor={args.zero_advective_flux_from_dry_donor}, "
          f"no_leak_rim_gate={args.no_leak_rim_gate}, scale={args.scale:g}: steps={steps}, "
          f"final_total_mass={total_mass:.12e} kg/m, "
          f"minimum_accepted_h={min_h:.12e} m, "
          f"Euler_retries={euler_retries}, "
          f"initial_M_uphill={initial_uphill_mass:.12e} kg/m, "
          f"max_uphill_mass_increase={max_uphill:.12e} kg/m")
    if local_blend:
        print("  initial alpha: "
              f"min={np.min(alpha0):.3g}, median={np.median(alpha0):.3g}, "
              f"max={np.max(alpha0):.3g}; "
              f"intermediate cells={np.sum((alpha0 > 1.e-12) & (alpha0 < 1.-1.e-12))}, "
              f"eta fallbacks={np.sum(infeasible0)}")
    for time, mass, advective, dissipative, left_adv, left_diss in reports:
        print(f"  t={time:5.2f} s  M_uphill={mass:.12e} kg/m  "
              f"face inflow (adv,diff)=({advective:+.8e}, "
              f"{dissipative:+.8e}) kg/m  "
              f"left-boundary contribution=({left_adv:+.8e}, "
              f"{left_diss:+.8e}) kg/m")


if __name__ == "__main__":
    main()
