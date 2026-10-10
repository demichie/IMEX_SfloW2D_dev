"""Characterize production IMEX stages, CFL policies and temporal order.

Reference tableaux are independent of private production storage. Spatial RHS
calls use the unchanged, checksum-pinned cores from TEST_SPATIAL_OPERATOR.
The thermal-mode and quadratic-drag order oracles are independent analytic
or method-of-lines calculations, not production fine-step trajectories.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "TEST_SPATIAL_OPERATOR"))
import compare_reference as spatial

CP = 4180.0
RHO = 1000.0
EPS = np.finfo(float).eps
LEVELS = (8, 16, 32, 64, 128)


def tableaux(n):
    """Return effective strictly lower explicit and lower implicit production tableaux."""
    ae = np.zeros((n, n)); ai = np.zeros((n, n))
    be = np.zeros(n); bi = np.zeros(n)
    if n == 1:
        # The stored explicit diagonal is unused by stage assembly.
        be[0] = 1.0
    elif n == 2:
        ae[1, 0] = 1.0; be[0] = 1.0; ai[1, 1] = 1.0; bi[1] = 1.0
    elif n == 3:
        ae[1, 0] = 0.5; ae[2, :2] = 0.5; be[:] = 1/3
        ai[0, 0] = ai[1, 1] = 0.25; ai[2, :] = 1/3; bi[:] = 1/3
    elif n == 4:
        ae[1, 0] = 0.5; ae[2, 0] = 1/3; ae[3, 1] = 1.0; be[1] = 1.0
        ai[1, 1] = 0.5; ai[2, 2] = 1/3; ai[3, 2:] = (0.75, 0.25); bi[2:] = (0.75, 0.25)
    else:
        raise AssertionError("unsupported tableau")
    return ae, ai, be, bi


def check_order_conditions():
    """Check first/second order and both mixed colored-tree conditions independently."""
    for n in range(1, 5):
        ae, ai, be, bi = tableaux(n)
        ce, ci = ae.sum(axis=1), ai.sum(axis=1)
        assert abs(be.sum()-1) < 8*EPS
        if n == 1:
            assert bi.sum() == 0  # No implicit source participates in this option.
        else:
            assert abs(bi.sum()-1) < 8*EPS
        if n >= 3:
            for weight in (be, bi):
                for abscissa in (ce, ci):
                    assert abs(weight@abscissa-0.5) < 8*EPS
        elif n == 2:
            assert be@ce != 0.5 and bi@ci != 0.5
    print("PASS: independent first/second and mixed order conditions; N_RK is stage count", flush=True)


def write_fixture(path, vertices, production, slope, curvature, limiter, dt, dx, dy):
    ny, nx, _ = production.shape
    with path.open("w") as stream:
        stream.write(f"{nx} {ny} {int(slope)} {int(curvature)} {limiter}\n{dx:.17g} {dy:.17g}\n{dt:.17g}\n")
        np.savetxt(stream, vertices.reshape(-1, 1), fmt="%.17g")
        np.savetxt(stream, production.reshape(-1, 1), fmt="%.17g")


def read_payload(path, nx, ny, n, steps):
    raw = np.fromfile(path, dtype=np.float64)
    fields = {}; cursor = 0

    def take(key, shape):
        nonlocal cursor
        count = int(np.prod(shape))
        fields[key] = raw[cursor:cursor+count].reshape(shape)
        cursor += count

    take("q0", (ny, nx, 5)); take("q", (ny, nx, 5)); take("qp", (ny, nx, 7))
    take("active", (ny, nx))
    take("known", (n, ny, nx, 5)); take("solved", (n, ny, nx, 5))
    take("raw_final", (ny, nx, 5)); take("implicit_statuses", (n, ny, nx))
    take("statistics", (steps, 15))
    if cursor != len(raw) or not np.all(np.isfinite(raw)):
        raise AssertionError("invalid IMEX payload length or nonfinite data")
    return fields


def launch(executable, directory, vertices, production, n, steps, dt, drag=0,
           slope=False, curvature=False, limiter=3, dx=0.75, dy=1.25, observe=1):
    directory.mkdir(parents=True)
    write_fixture(directory/"fixture.inp", vertices, production, slope, curvature, limiter, dt, dx, dy)
    results = []
    for threads in (1, 4):
        work = directory/f"threads-{threads}"; work.mkdir()
        (work/"fixture.inp").write_bytes((directory/"fixture.inp").read_bytes())
        env = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE"}
        with (work/"solver.log").open("w") as stream:
            run = subprocess.run([str(executable), str(threads), str(n), str(steps), str(drag), str(observe)],
                                 cwd=work, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=180)
        if run.returncode:
            raise AssertionError(f"{work}: exit {run.returncode}\n"+(work/"solver.log").read_text())
        ny, nx, _ = production.shape
        results.append(read_payload(work/"result.bin", nx, ny, n, steps))
    first, second = directory/"threads-1/result.bin", directory/"threads-4/result.bin"
    if first.read_bytes() != second.read_bytes():
        raise AssertionError(f"{directory}: actual 1/4-thread outputs differ")
    evidence = {"input_sha256": hashlib.sha256((directory/"fixture.inp").read_bytes()).hexdigest(),
                "output_sha256": hashlib.sha256(first.read_bytes()).hexdigest(), "actual_threads": [1, 4]}
    return results[0], evidence


def close(label, actual, expected, scale=None, mask=None):
    actual, expected = np.broadcast_arrays(actual, expected)
    if mask is not None:
        actual, expected = actual[mask], expected[mask]
    scale = max(1.0, float(np.max(np.abs(expected)))) if scale is None else max(1.0, scale)
    error = float(np.max(np.abs(actual-expected)))/scale
    if error > 32768*EPS:
        raise AssertionError(f"{label}: scaled error {error:.17g} > {32768*EPS:.17g}")
    return error


def admissibility(fields):
    """Inspect raw conservative states, not clipped/desingularized primitive values."""
    stats = fields["statistics"]
    if np.any(stats[:, 2:5] < -32768*EPS):
        raise AssertionError("negative raw stage/final thickness")
    if np.any(stats[:, 6] != 0):
        raise AssertionError("failed implicit Newton solve")
    if np.any(stats[:, 12] < 273-32768*EPS) or np.any(stats[:, 13:15] < -32768*EPS):
        raise AssertionError("inadmissible raw temperature/component/carrier fraction")
    active = fields["active"] > 0
    for key in ("known", "solved", "raw_final"):
        states = fields[key][..., active, :] if fields[key].ndim == 4 else fields[key][active]
        mass = states[..., 0]; wet = mass > 0
        if np.any(states[..., 4] < -32768*EPS) or np.any(mass-states[..., 4] < -32768*EPS):
            raise AssertionError("inadmissible raw component/carrier mass")
        if np.any(states[..., 3][wet]/(CP*mass[wet]) < 273-32768*EPS):
            raise AssertionError("inadmissible raw temperature")
    return {"minimum_known_h": float(stats[:, 2].min()),
            "minimum_solved_h": float(stats[:, 3].min()), "minimum_raw_final_h": float(stats[:, 4].min()),
            "maximum_final_repair": float(stats[:, 5].max()), "failed_local_solves": int(stats[:, 6].sum()),
            "converged_local_solves": int(stats[:, 7].sum()), "minimum_raw_temperature": float(stats[:, 12].min()),
            "minimum_solid_fraction": float(stats[:, 13].min()), "minimum_carrier_fraction": float(stats[:, 14].min())}


def lift(q):
    production = np.zeros(q.shape[:-1]+(5,))
    production[..., :3] = q[..., :3]; production[..., 3] = spatial.THERMAL_FACTOR*q[..., 0]
    production[..., 4] = q[..., 3]
    return production


def explicit_stages(q, bed, params, n, dt):
    ae, _, be, _ = tableaux(n)
    terms = []; states = []; stage_cfl = []
    for i in range(n):
        state = q.copy()
        for j in range(i):
            state += dt*ae[i, j]*terms[j]
        states.append(lift(state))
        cache = spatial.core.prepare_stage_cache2d(state, bed, 0.0, params, dx=spatial.DX, dy=spatial.DY)
        terms.append(spatial.core.rhs2d(state, bed, spatial.DX, spatial.DY, 0.0, params, stage_cache=cache))
        harmonic, bounds = spatial.core.max_dt2d(state, bed, spatial.DX, spatial.DY, 0.0, params, stage_cache=cache)
        stage_cfl.append({"stage": i+1, "dt_x": bounds["dt_x"], "dt_y": bounds["dt_y"],
                          "dt_harmonic": harmonic, "dt_over_directional_minimum": dt/min(bounds["dt_x"], bounds["dt_y"]),
                          "combined_Courant": dt*params.cfl/harmonic})
    final = q.copy()
    for weight, rhs in zip(be, terms):
        final += dt*weight*rhs
    return np.array(states), lift(final), stage_cfl


def validate_stages(fields, expected, final, dt, dt_bound, lake=False):
    active = fields["active"] > 0
    for state in expected:
        if np.any(~active):
            close("reference outside active workset", state[~active], 0.0)
    errors = {}
    # Do not let the much larger thermal-energy units hide a momentum error.
    for component, label in enumerate(("mass", "mx", "my", "thermal", "solid")):
        for key in ("known", "solved"):
            errors[key+"_"+label] = close(key+"_"+label, fields[key][..., component], expected[..., component],
                                        mask=np.broadcast_to(active, expected.shape[:-1]))
        errors["raw_final_"+label] = close("raw_final_"+label, fields["raw_final"][..., component],
                                          final[..., component], mask=active)
        errors["final_"+label] = close("final_"+label, fields["q"][..., component], final[..., component])
    errors["dt_used"] = close("dt used", fields["statistics"][0, 0], dt)
    errors["dt_bound"] = close("dt bound", fields["statistics"][0, 1], dt_bound)
    errors["mass_budget"] = close("global mass budget", fields["statistics"][:, 8:10], 0.0,
                                 scale=float(np.sum(np.abs(fields["q0"][..., 0]))))
    if lake:
        errors["lake_unchanged"] = close("lake at rest advanced", fields["q"], fields["q0"])
    return errors


def negative_controls(fields, expected, final, dt, bound):
    for key in ("known", "solved", "raw_final"):
        bad = {k: v.copy() for k, v in fields.items()}
        index = (0, 8, 9, 1) if bad[key].ndim == 4 else (8, 9, 1)
        bad[key][index] += 1.0
        try:
            validate_stages(bad, expected, final, dt, bound)
        except AssertionError as error:
            assert key in str(error)
        else:
            raise AssertionError(f"undetected {key} corruption")
    for column, amount, label in ((2, -1.0, "thickness"), (6, 1.0, "Newton"),
                                  (12, -300.0, "temperature"), (14, -2.0, "fraction")):
        bad = {k: v.copy() for k, v in fields.items()}; bad["statistics"][0, column] += amount
        try:
            admissibility(bad)
        except AssertionError as error:
            assert label in str(error)
        else:
            raise AssertionError(f"undetected {label} corruption")
    print("PASS: seven raw-stage comparator/admissibility negative controls", flush=True)


def thermal_rhs(temperature, u, v, dx=0.75, dy=1.25):
    """Independent first-order CU thermal operator, with Neumann ghost temperatures."""
    c = np.sqrt(9.81*(1-1.2/RHO))
    left = np.concatenate((temperature[:, :1], temperature), axis=1)
    right = np.concatenate((temperature, temperature[:, -1:]), axis=1)
    bottom = np.concatenate((temperature[:1, :], temperature), axis=0)
    top = np.concatenate((temperature, temperature[-1:, :]), axis=0)
    fx = 0.5*(u+c)*left + 0.5*(u-c)*right
    fy = 0.5*(v+c)*bottom + 0.5*(v-c)*top
    return -(fx[:, 1:]-fx[:, :-1])/dx-(fy[1:]-fy[:-1])/dy


def thermal_oracle(initial, final_time, drag, steps):
    """RK4 of an independent linear thermal operator with analytic uniform drag velocity."""
    state = initial.copy(); dt = final_time/steps

    def rhs(t, state):
        speed = 1/(1+t) if drag else 1.0
        return thermal_rhs(state, 0.6*speed, 0.8*speed)

    for i in range(steps):
        t = i*dt
        k1 = rhs(t, state); k2 = rhs(t+dt/2, state+dt*k1/2)
        k3 = rhs(t+dt/2, state+dt*k2/2); k4 = rhs(t+dt, state+dt*k3)
        state += dt*(k1+2*k2+2*k3+k4)/6
    return state


def manufactured(kind):
    ny, nx = 7, 8
    x, y = np.meshgrid(np.arange(nx)+0.5, np.arange(ny)+0.5)
    mode = np.cos(np.pi*x/nx)*np.cos(np.pi*y/ny)
    temperature = 300+15*mode if kind != "drag" else np.full((ny, nx), 300.0)
    q = np.zeros((ny, nx, 5)); q[..., 0] = RHO; q[..., 3] = RHO*CP*temperature
    q[..., 1] = 600; q[..., 2] = 800
    return np.zeros((ny+1, nx+1)), q, temperature, mode


def implicit_stage_checks(fields, n, dt):
    """Check real nonzero drag solves against their quadratic equation and root.

    Production accepts either its residual test or a small normalized Newton
    correction. Check both alternatives rather than inventing a residual-only
    contract or claiming machine-precision implicit solutions.
    """
    _, ai, _, bi = tableaux(n)
    known, solved = fields["known"], fields["solved"]
    increments = []; largest_residual = 0.0; largest_root_error = 0.0
    correction_only = 0; largest_relative_correction = 0.0
    start = known[0, ..., 1:3]
    for i in range(n):
        expected_known = start.copy()
        for j in range(i):
            expected_known += dt*ai[i, j]*increments[j]
        close("implicit known momentum assembly", known[i, ..., 1:3], expected_known)
        mass = solved[i, ..., 0]
        # h=1 in this manufactured fixture: m'=-|m|m/rho.
        actual = solved[i, ..., 1:3]; old = known[i, ..., 1:3]
        force = -np.linalg.norm(actual, axis=-1)[..., None]*actual/RHO
        diagonal = ai[i, i]
        if diagonal == 0:
            close("zero diagonal implicit stage", actual, old)
            increments.append(np.zeros_like(actual))
            continue
        residue = actual-old-dt*diagonal*force
        initial_force = -np.linalg.norm(old, axis=-1)[..., None]*old/RHO
        tolerance = 1e-5+1e-5*np.abs(dt*diagonal*initial_force)
        residual_met = np.all(np.abs(residue) <= tolerance+32768*EPS*np.maximum(1, np.abs(old)), axis=-1)
        # For uniform collinear drag the Newton correction is radial. The
        # production normalized variables use q_org=max(abs(initial_guess),1e-3).
        correction = residue/(1+2*dt*diagonal*np.linalg.norm(actual, axis=-1)[..., None]/RHO)
        relative = np.max(np.abs(correction)/np.maximum(1e-3, np.abs(old)), axis=-1)
        if np.any(~residual_met & (relative > 1e-5+32768*EPS)):
            raise AssertionError("quadratic drag fails both production Newton convergence criteria")
        correction_only += int(np.sum(~residual_met))
        largest_relative_correction = max(largest_relative_correction, float(relative.max()))
        old_norm = np.linalg.norm(old, axis=-1)
        root = 2*old_norm/(1+np.sqrt(1+4*dt*diagonal*old_norm/RHO))
        exact = old*(root/old_norm)[..., None]
        root_limit = np.maximum(2*tolerance, 2e-5*np.maximum(1e-3, np.abs(old)))
        if np.any(np.abs(actual-exact) > root_limit):
            raise AssertionError("quadratic drag stage differs from its independent algebraic root")
        largest_residual = max(largest_residual, float(np.max(np.abs(residue))))
        largest_root_error = max(largest_root_error, float(np.max(np.abs(actual-exact))))
        increments.append((actual-old)/(dt*diagonal))
    final = start.copy()
    for weight, increment in zip(bi, increments):
        final += dt*weight*increment
    close("implicit final momentum assembly", fields["raw_final"][..., 1:3], final)
    return {"maximum_last_step_implicit_residual": largest_residual,
            "maximum_last_step_quadratic_root_error": largest_root_error,
            "last_step_cells_using_correction_only_criterion": correction_only,
            "maximum_last_step_relative_Newton_correction": largest_relative_correction,
            "production_Newton_tol_abs": 1e-5, "production_Newton_tol_rel": 1e-5}


def order_suite(args):
    final_time = 0.32
    for kind in ("thermal", "drag", "mixed"):
        vertices, q0, initial_t, mode = manufactured(kind)
        if kind == "thermal":
            # Exact semidiscrete matrix exponential, including the Neumann
            # boundary rows. A resting thermal mode would trigger the intended
            # HP equilibrium scalar-flux suppression and is not an order test.
            count = initial_t.size
            basis = np.eye(count).reshape(count, *initial_t.shape)
            matrix = np.column_stack([thermal_rhs(e, 0.6, 0.8).ravel() for e in basis])
            eigenvalues, vectors = np.linalg.eig(matrix)
            expected_t = vectors @ (np.exp(final_time*eigenvalues)*np.linalg.solve(vectors, initial_t.ravel()))
            if np.max(np.abs(expected_t.imag)) > 32768*EPS*np.max(np.abs(expected_t.real)):
                raise AssertionError("complex contamination in exact thermal oracle")
            expected_t = expected_t.real.reshape(initial_t.shape)
            close("independent thermal oracle versus exact MOL exponential",
                  thermal_oracle(initial_t, final_time, False, 2048), expected_t)
        elif kind == "mixed":
            expected_t = thermal_oracle(initial_t, final_time, True, 4096)
            close("mixed oracle refinement", thermal_oracle(initial_t, final_time, True, 2048), expected_t)
        else:
            expected_t = initial_t
        for n in range(2, 5):
            errors = []; fingerprints = []; minima = []; stage_checks = []
            for steps in LEVELS:
                directory = Path(f"order-{kind}-RK{n}-steps{steps}")
                fields, hashes = launch(args.executable, directory, vertices, q0, n, steps,
                                        final_time/steps, drag=int(kind != "thermal"), limiter=0)
                minima.append(admissibility(fields)); fingerprints.append(hashes)
                if np.any(fields["statistics"][:, 0] > fields["statistics"][:, 1]*(1+32768*EPS)):
                    raise AssertionError("order fixture exceeds the retained production CFL bound")
                close("manufactured mass remains constant", fields["q"][..., 0], RHO)
                close("manufactured solid remains zero", fields["q"][..., 4], 0.0)
                close("no raw-stage clipping in order fixture", fields["statistics"][:, 5], 0.0)
                if kind != "thermal":
                    stage_checks.append(implicit_stage_checks(fields, n, final_time/steps))
                    expected_momentum = q0[..., 1:3]/(1+final_time)
                    error = np.max(np.abs(fields["q"][..., 1:3]-expected_momentum))/RHO
                    if kind == "mixed":
                        error = max(error, float(np.max(np.abs(fields["qp"][..., 3]-expected_t))))
                else:
                    error = np.max(np.abs(fields["qp"][..., 3]-expected_t))
                errors.append(float(error))
            rates = np.log2(np.array(errors[:-1])/errors[1:])
            target = 1 if n == 2 or n == 1 else 2
            if not np.all((rates[-2:] > target-0.25) & (rates[-2:] < target+0.25)):
                raise AssertionError(f"temporal order {kind} RK{n}: errors={errors}, rates={rates}")
            record = {"profile": args.profile, "kind": kind, "n_RK": n, "expected_order": target,
                      "final_time": final_time, "step_counts": LEVELS, "errors": errors,
                      "observed_orders": rates.tolist(), "fingerprints": fingerprints, "admissibility": minima,
                      "implicit_stage_checks": stage_checks,
                      "oracle": "exact constant-coefficient thermal MOL exponential" if kind == "thermal" else
                      "analytic quadratic drag" if kind == "drag" else "independent thermal MOL RK4 plus analytic drag"}
            print("IMEX_ORDER_EVIDENCE "+json.dumps(record, sort_keys=True, allow_nan=False), flush=True)
    # N_RK=1 has zero implicit weights: it must not be claimed to integrate drag.
    vertices, q0, _, _ = manufactured("drag")
    fields, hashes = launch(args.executable, Path("RK1-implicit-disabled"), vertices, q0, 1, 8, 0.05, drag=1, limiter=0)
    close("N_RK=1 does not integrate implicit friction", fields["q"], q0)
    if np.any(fields["implicit_statuses"] != 0):
        raise AssertionError("unexpected N_RK=1 implicit solve")
    print("PASS: N_RK=1 implicit contribution is disabled by its retained tableau", flush=True)
    # Reproduce the separate explicit defect without declaring Euler agreement.
    vertices, q0, initial_t, mode = manufactured("thermal")
    fields, hashes = launch(args.executable, Path("RK1-explicit-noop"), vertices, q0, 1, 8, 0.05, limiter=0)
    if fields["q"].tobytes() != q0.tobytes():
        raise AssertionError("N_RK=1 known no-op changed; update the explicit defect gate")
    print("IMEX_KNOWN_DEFECT "+json.dumps({"profile": args.profile, "n_RK": 1,
          "defect": "stored explicit diagonal incorrectly selects stiffly-accurate assembly; explicit update is skipped",
          "intended_forward_Euler_matches": False, **hashes}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path); parser.add_argument("profile")
    parser.add_argument("--section", choices=("stages", "order"))
    parser.add_argument("--require-design-order", action="store_true",
                        help="fail on characterized numerical defects, not just unexpected regressions")
    args = parser.parse_args(); args.executable = args.executable.resolve()
    check_order_conditions()
    if args.section != "order":
        cases = ("lake_rough", "dynamic_smooth", "dynamic_drainage", "dynamic_compact")
        for name in cases:
            vertices, bed, q = spatial.fixture(name)
            for slope, curvature in ((False, False), (True, True)):
                params = spatial.core.SolverParams(rho_a=1.2, eps_sing=1e-8,
                                                   slope_correction=slope, curvature_term=curvature)
                harmonic, diagnostic = spatial.core.max_dt2d(q, bed, spatial.DX, spatial.DY, 0.0, params)
                bound = min(diagnostic["dt_x"], diagnostic["dt_y"])
                for n in range(1, 5):
                    for policy, requested, dt in (("retained", 0.0, bound), ("harmonic", harmonic, harmonic)):
                        expected, final, stage_cfl = explicit_stages(q, bed, params, n, dt)
                        tag = f"{name}-G{int(slope)}-C{int(curvature)}-RK{n}-{policy}"
                        fields, hashes = launch(args.executable, Path(tag), vertices, lift(q), n, 1, requested,
                                                slope=slope, curvature=curvature, limiter=params.limiter)
                        intended_final = final
                        if n == 1:
                            # Keep the reference Euler result as evidence; production
                            # currently returns Q^n due to its diagonal/stiff-accuracy test.
                            final = lift(q)
                            if fields["q"].tobytes() != lift(q).tobytes():
                                raise AssertionError("N_RK=1 no-op defect changed; re-evaluate its gate")
                        errors = validate_stages(fields, expected, final, dt, bound, name.startswith("lake"))
                        record = {"profile": args.profile, "case": tag, "n_RK": n, "policy": policy,
                                  "dt_used": dt, "dt_directional_minimum": bound, "dt_harmonic": harmonic,
                                  "active_cells": int(fields["active"].sum()), "maximum_scaled_errors": errors,
                                  "stage_CFL_diagnostics": stage_cfl,
                                  "admissibility": admissibility(fields), **hashes}
                        record["intended_final_component_scaled_errors"] = {
                            label: float(np.max(np.abs(fields["raw_final"][..., i]-intended_final[..., i]))) /
                            max(1.0, float(np.max(np.abs(intended_final[..., i]))))
                            for i, label in enumerate(("mass", "mx", "my", "thermal", "solid"))}
                        record["known_RK1_noop"] = n == 1
                        print("IMEX_STAGE_EVIDENCE "+json.dumps(record, sort_keys=True, allow_nan=False), flush=True)
                        if name == "dynamic_smooth" and slope and n == 3 and policy == "retained":
                            negative_controls(fields, expected, final, dt, bound)
                        if name in ("lake_rough", "dynamic_smooth") and slope and policy == "retained":
                            without, _ = launch(args.executable, Path(tag+"-observer-off"), vertices, lift(q), n, 1,
                                                requested, slope=slope, curvature=curvature, observe=0)
                            for key in ("q", "qp"):
                                if without[key].tobytes() != fields[key].tobytes():
                                    raise AssertionError("optional observer changed the production state")
        print("PASS: production raw stages and assembly, both CFL policies, actual 1/4 threads", flush=True)
    if args.section != "stages":
        order_suite(args)
    if args.require_design_order:
        raise AssertionError("N_RK=1 skips its intended explicit update; design-order acceptance remains open")


if __name__ == "__main__":
    main()
