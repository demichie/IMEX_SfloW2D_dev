"""Validate supported IMEX stages, CFL policies, temporal order and input rejection.

Reference tableaux are independent of private production storage. Spatial RHS
calls use the checksum-pinned independent KT reference from TEST_SPATIAL_OPERATOR.
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
TEST_STAGES = (2,)
TABLEAU_STAGES = (2, 3, 4)
EPS = np.finfo(float).eps
LEVELS = (8, 16, 32, 64, 128)


def tableaux(n):
    """Return effective strictly lower explicit and lower implicit production tableaux."""
    ae = np.zeros((n, n)); ai = np.zeros((n, n))
    be = np.zeros(n); bi = np.zeros(n)
    if n == 2:
        ae[1, 0] = 1.0; be[0] = 1.0; ai[1, 1] = 1.0; bi[1] = 1.0
    elif n == 3:
        ae[1, 0] = 0.5; ae[2, :2] = 0.5; be[:] = 1/3
        ai[0, 0] = ai[1, 1] = 0.25; ai[2, :] = 1/3; bi[:] = 1/3
    elif n == 4:
        # IMEX-SSP3(4,3,3), independently restated from the published
        # tableau used by the production n_RK=4 implementation.
        gamma = 0.24169426078821
        delta = 0.06042356519705
        ae[2, 1] = 1.0
        ae[3, 1:3] = (0.25, 0.25)
        be[1:] = (1/6, 1/6, 2/3)
        ai[0, 0] = gamma
        ai[1, 0:2] = (-gamma, gamma)
        ai[2, 1:3] = (1.0-gamma, gamma)
        ai[3, :] = (delta, 0.12915286960590,
                    0.5-delta-0.12915286960590-gamma, gamma)
        bi[:] = be
    else:
        raise AssertionError("unsupported tableau")
    return ae, ai, be, bi


def check_order_conditions():
    """Check first/second order and both mixed colored-tree conditions independently."""
    for n in TABLEAU_STAGES:
        ae, ai, be, bi = tableaux(n)
        ce, ci = ae.sum(axis=1), ai.sum(axis=1)
        assert abs(be.sum()-1) < 8*EPS
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
           slope=False, curvature=False, limiter=3, dx=0.75, dy=1.25, observe=1,
           analytic_mode=0, source_test_model=0, mixed_group_mode=0):
    directory.mkdir(parents=True)
    write_fixture(directory/"fixture.inp", vertices, production, slope, curvature, limiter, dt, dx, dy)
    results = []
    for threads in (1, 4):
        work = directory/f"threads-{threads}"; work.mkdir()
        (work/"fixture.inp").write_bytes((directory/"fixture.inp").read_bytes())
        env = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE"}
        with (work/"solver.log").open("w") as stream:
            command = [str(executable), str(threads), str(n), str(steps), str(drag), str(observe)]
            if analytic_mode:
                command.extend(("0", "0", "0", str(analytic_mode)))
            elif mixed_group_mode:
                command.extend(("0", "0", str(source_test_model), "0", str(mixed_group_mode)))
            run = subprocess.run(command, cwd=work, env=env, stdout=stream,
                                 stderr=subprocess.STDOUT, timeout=180)
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


def explicit_stages(q, vertices, params, n, dt):
    ae, _, be, _ = tableaux(n)
    terms = []; states = []; stage_cfl = []
    for i in range(n):
        state = q.copy()
        for j in range(i):
            state += dt*ae[i, j]*terms[j]
        states.append(lift(state))
        term, diagnostics = spatial.kt.rhs(
            state, vertices, spatial.DX, spatial.DY, params, return_diagnostics=True)
        terms.append(term)
        dt_x = params.cfl*spatial.DX/float(np.max(np.maximum(
            -diagnostics["x"]["a_minus"], diagnostics["x"]["a_plus"])))
        dt_y = params.cfl*spatial.DY/float(np.max(np.maximum(
            -diagnostics["y"]["a_minus"], diagnostics["y"]["a_plus"])))
        harmonic = 1.0/(1.0/dt_x+1.0/dt_y)
        stage_cfl.append({"stage": i+1, "dt_x": dt_x, "dt_y": dt_y,
                          "dt_harmonic": harmonic, "dt_over_directional_minimum": dt/min(dt_x, dt_y),
                          "combined_Courant": dt/harmonic})
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
    ny, nx = 16, 16
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
    _, ai, _, _ = tableaux(n)
    known, solved = fields["known"], fields["solved"]
    largest_residual = 0.0; largest_root_error = 0.0
    correction_only = 0; largest_relative_correction = 0.0
    for i in range(n):
        actual = solved[i, ..., 1:3]; old = known[i, ..., 1:3]
        depth = known[i, ..., 0]/RHO
        depth_squared = depth**2
        force = -np.linalg.norm(actual, axis=-1)[..., None]*actual / (
            RHO*depth_squared[..., None])
        diagonal = ai[i, i]
        if diagonal == 0:
            close("zero diagonal implicit stage", actual, old)
            continue
        residue = actual-old-dt*diagonal*force
        initial_force = -np.linalg.norm(old, axis=-1)[..., None]*old / (
            RHO*depth_squared[..., None])
        tolerance = 1e-5+1e-5*np.abs(dt*diagonal*initial_force)
        residual_met = np.all(np.abs(residue) <= tolerance+32768*EPS*np.maximum(1, np.abs(old)), axis=-1)
        # For uniform collinear drag the Newton correction is radial. The
        # production normalized variables use q_org=max(abs(initial_guess),1e-3).
        correction = residue/(1+2*dt*diagonal*np.linalg.norm(actual, axis=-1)[..., None]/
                              (RHO*depth_squared[..., None]))
        relative = np.max(np.abs(correction)/np.maximum(1e-3, np.abs(old)), axis=-1)
        if np.any(~residual_met & (relative > 1e-5+32768*EPS)):
            raise AssertionError("quadratic drag fails both production Newton convergence criteria")
        correction_only += int(np.sum(~residual_met))
        largest_relative_correction = max(largest_relative_correction, float(relative.max()))
        old_norm = np.linalg.norm(old, axis=-1)
        root = 2*old_norm/(1+np.sqrt(
            1+4*dt*diagonal*old_norm/(RHO*depth_squared)))
        exact = old*(root/old_norm)[..., None]
        root_limit = np.maximum(2*tolerance, 2e-5*np.maximum(1e-3, np.abs(old)))
        if np.any(np.abs(actual-exact) > root_limit):
            raise AssertionError("quadratic drag stage differs from its independent algebraic root")
        largest_residual = max(largest_residual, float(np.max(np.abs(residue))))
        largest_root_error = max(largest_root_error, float(np.max(np.abs(actual-exact))))
    return {"maximum_last_step_implicit_residual": largest_residual,
            "maximum_last_step_quadratic_root_error": largest_root_error,
            "last_step_cells_using_correction_only_criterion": correction_only,
            "maximum_last_step_relative_Newton_correction": largest_relative_correction,
            "production_Newton_tol_abs": 1e-5, "production_Newton_tol_rel": 1e-5}


def analytic_quadratic_stage_checks(fields, n, dt, vertices):
    """Check ANALYTIC stages against independent frozen-coefficient drag flows."""
    ae, ai, be, bi = tableaux(n)
    start = fields["known"][0, ..., 1:3]
    depth = fields["q0"][..., 0]/RHO
    params = spatial.kt.Parameters(rho_a=1.2, eps_sing=1e-8, limiter=0)
    interior = np.zeros(start.shape[:2], dtype=bool)
    interior[4:-4, 4:-4] = True
    explicit_rates = []
    increments = []
    for i in range(n):
        expected_known = start.copy()
        for j in range(i):
            expected_known += dt*ae[i, j]*explicit_rates[j][..., 1:3]
            expected_known += dt*ai[i, j]*increments[j]
        close("analytic known momentum assembly", fields["known"][i, ..., 1:3], expected_known,
              mask=np.broadcast_to(interior[..., None], expected_known.shape))
        stage = fields["solved"][i]
        physical = np.stack((stage[..., 0], stage[..., 1], stage[..., 2], stage[..., 4]), axis=-1)
        explicit_rates.append(spatial.kt.rhs(physical, vertices, 0.75, 1.25, params))
        diagonal = ai[i, i]
        old = fields["known"][i, ..., 1:3]
        actual = fields["solved"][i, ..., 1:3]
        if diagonal == 0:
            close("zero diagonal analytic stage", actual, old,
                  mask=np.broadcast_to(interior[..., None], old.shape))
            increments.append(np.zeros_like(actual))
            continue
        magnitude = np.linalg.norm(old, axis=-1)
        exact = old/(1+dt*diagonal*magnitude/(RHO*depth**2))[..., None]
        close("exact-flow analytic stage", actual, exact,
              mask=np.broadcast_to(interior[..., None], exact.shape))
        increments.append((actual-old)/(dt*diagonal))
    final = start.copy()
    for weight, increment in zip(bi, increments):
        final += dt*weight*increment
    for weight, rate in zip(be, explicit_rates):
        final += dt*weight*rate[..., 1:3]
    close("analytic final momentum assembly", fields["raw_final"][..., 1:3], final,
          mask=np.broadcast_to(interior[..., None], final.shape))
    close("analytic repaired momentum assembly", fields["q"][..., 1:3], final,
          mask=np.broadcast_to(interior[..., None], final.shape))
    exact_final = start/(1+dt*np.linalg.norm(start, axis=-1)/(RHO*depth**2))[..., None]
    return float(np.max(np.abs(fields["q"][..., 1:3][interior] -
                               exact_final[interior])))/RHO


def analytic_tableau_suite(args):
    """Exercise all supported tableaux through the production ANALYTIC dispatcher."""
    ny, nx = 17, 18
    vertices = np.zeros((ny+1, nx+1))
    q0 = np.zeros((ny, nx, 5))
    q0[..., 0] = RHO
    q0[..., 1] = 600.0
    q0[..., 2] = 800.0
    q0[..., 3] = RHO*CP*300.0
    bed = np.zeros(q0.shape[:2])
    dt_values = (0.08, 0.04, 0.02)
    for n in TABLEAU_STAGES:
        errors = []
        for dt in dt_values:
            directory = Path(f"analytic-{args.profile}-{os.getpid()}-quadratic-RK{n}-dt{dt:g}")
            fields, hashes = launch(args.executable, directory, vertices, q0, n, 1, dt,
                                    drag=1, limiter=0, analytic_mode=1)
            interior = np.zeros(fields["active"].shape, dtype=bool)
            interior[4:-4, 4:-4] = True
            close("analytic fixture mass", fields["q"][..., 0], RHO, mask=interior)
            close("analytic fixture solid", fields["q"][..., 4], 0.0, mask=interior)
            if np.any(fields["statistics"][:, 5] != 0):
                raise AssertionError("analytic tableau fixture unexpectedly clipped a raw stage")
            error = analytic_quadratic_stage_checks(fields, n, dt, vertices)
            errors.append(error)
            print("IMEX_ANALYTIC_TABLEAU_EVIDENCE "+json.dumps({
                "profile": args.profile, "n_RK": n, "dt": dt,
                "exact_quadratic_endpoint_error": error, **hashes
            }, sort_keys=True, allow_nan=False), flush=True)
        if errors[0] <= 1024*EPS:
            if max(errors) > 1024*EPS:
                raise AssertionError(f"exact-flow endpoint lost exactness under refinement: RK{n}, {errors}")
        elif not errors[-1] < errors[0]:
            raise AssertionError(f"exact-flow tableau does not converge under dt refinement: RK{n}, {errors}")
    print("PASS: exact-flow stage trajectories and final assembly for N_RK=2/3/4", flush=True)


def close_approx(label, actual, expected, relative_tolerance, mask):
    actual, expected = np.broadcast_arrays(actual, expected)
    actual, expected = actual[mask], expected[mask]
    scale = max(1.0, float(np.max(np.abs(expected))))
    error = float(np.max(np.abs(actual-expected)))/scale
    if error > relative_tolerance:
        raise AssertionError(f"{label}: scaled error {error:.17g} > {relative_tolerance:.17g}")


def newton_lava_temperature(energy, mass, interval):
    """Independent bisection of the scalar DIRK cooling residual."""
    ambient = 300.0
    heat_capacity = CP*mass
    initial_temperature = energy/heat_capacity
    lower = np.full_like(initial_temperature, ambient)
    upper = initial_temperature.copy()
    for _ in range(80):
        midpoint = 0.5*(lower+upper)
        source = (-3.0e-8*(midpoint**4-ambient**4) -
                  1.0e4*(midpoint-ambient))
        residual = heat_capacity*(midpoint-initial_temperature)-interval*source
        lower = np.where(residual < 0.0, midpoint, lower)
        upper = np.where(residual >= 0.0, midpoint, upper)
    return 0.5*(lower+upper)


def exact_lava_momentum(momentum, mass, temperature, interval):
    """Independent frozen-temperature affine-drag velocity flow."""
    depth = mass/RHO
    speed = np.linalg.norm(momentum, axis=-1)/mass
    linear_drag = 3.0*0.02/depth**2*np.exp(-0.3*(temperature-1000.0))
    yield_drag = 20.0/(RHO*depth)
    speed_after = np.maximum(
        0.0, (speed+yield_drag/linear_drag)*np.exp(-linear_drag*interval) -
        yield_drag/linear_drag)
    return momentum*(speed_after/np.maximum(speed, np.finfo(float).tiny))[..., None]


def mixed_group_stage_checks(fields, n, dt, group_order):
    """Validate mixed ANALYTIC/NEWTON stages and tableau assembly in the interior."""
    _, ai, _, bi = tableaux(n)
    start = fields["q0"][..., :4]
    interior = np.zeros(start.shape[:2], dtype=bool)
    interior[4:-4, 4:-4] = True
    rates = []
    if np.any(fields["statistics"][:, 6] != 0.0) or fields["statistics"][0, 7] <= 0.0:
        raise AssertionError("mixed-group trajectory had a failed or missing Newton solve")
    for i in range(n):
        expected_known = start.copy()
        for j in range(i):
            expected_known += dt*ai[i, j]*rates[j]
        close("mixed-group known tableau assembly", fields["known"][i, ..., :4],
              expected_known, mask=np.broadcast_to(interior[..., None], expected_known.shape))
        known = fields["known"][i, ..., :4]
        solved = fields["solved"][i, ..., :4]
        diagonal = ai[i, i]
        if diagonal == 0.0:
            close("mixed-group zero-diagonal state", solved, known,
                  mask=np.broadcast_to(interior[..., None], known.shape))
            rates.append(np.zeros_like(known))
            continue

        interval = dt*diagonal
        mass = known[..., 0]
        known_temperature = known[..., 3]/(CP*mass)
        expected_temperature = newton_lava_temperature(known[..., 3], mass, interval)
        actual_temperature = solved[..., 3]/(CP*mass)
        temperature_error = float(np.max(np.abs(
            actual_temperature[interior]-expected_temperature[interior])))
        if temperature_error > 2.0e-5:
            raise AssertionError(
                f"mixed-group Newton cooling temperature error {temperature_error:.17g}")
        velocity_temperature = known_temperature if group_order == 1 else expected_temperature
        expected_momentum = exact_lava_momentum(
            known[..., 1:3], mass, velocity_temperature, interval)
        close_approx("mixed-group exact-flow momentum", solved[..., 1:3],
                     expected_momentum, 2.0e-10,
                     np.broadcast_to(interior[..., None], expected_momentum.shape))
        rates.append((solved-known)/(dt*diagonal))

    expected_final = start.copy()
    for weight, rate in zip(bi, rates):
        expected_final += dt*weight*rate
    close("mixed-group raw tableau assembly", fields["raw_final"][..., :4],
          expected_final, mask=np.broadcast_to(interior[..., None], expected_final.shape))
    close("mixed-group repaired tableau assembly", fields["q"][..., :4],
          expected_final, mask=np.broadcast_to(interior[..., None], expected_final.shape))
    return fields["q"][..., 1:3][interior].copy()


def mixed_group_tableau_suite(args):
    """Exercise both mixed-method orderings through every supported tableau."""
    ny, nx = 17, 18
    vertices = np.zeros((ny+1, nx+1))
    q0 = np.zeros((ny, nx, 5))
    q0[..., 0] = RHO
    q0[..., 1] = 600.0
    q0[..., 2] = 800.0
    q0[..., 3] = RHO*CP*1000.0
    for n in TABLEAU_STAGES:
        endpoints = {}
        for group_order in (1, 2):
            directory = Path(
                f"mixed-{args.profile}-{os.getpid()}-RK{n}-order{group_order}")
            fields, hashes = launch(
                args.executable, directory, vertices, q0, n, 1, 0.5,
                limiter=0, dx=100.0, dy=100.0, source_test_model=3,
                mixed_group_mode=group_order)
            interior = np.zeros(fields["active"].shape, dtype=bool)
            interior[4:-4, 4:-4] = True
            close("mixed-group fixture mass", fields["q"][..., 0], RHO, mask=interior)
            close("mixed-group fixture solid", fields["q"][..., 4], 0.0, mask=interior)
            endpoints[group_order] = mixed_group_stage_checks(
                fields, n, fields["statistics"][0, 0], group_order)
            print("IMEX_MIXED_GROUP_EVIDENCE "+json.dumps({
                "profile": args.profile, "n_RK": n, "group_order": group_order,
                "dt": float(fields["statistics"][0, 0]), **hashes
            }, sort_keys=True, allow_nan=False), flush=True)
        order_difference = float(np.max(np.abs(endpoints[1]-endpoints[2])))/RHO
        if order_difference <= 1.0e-8:
            raise AssertionError(
                f"configured mixed-group order did not affect the frozen-flow trajectory: "
                f"RK{n}, difference={order_difference:.17g}")
        print("IMEX_MIXED_GROUP_ORDER_EVIDENCE "+json.dumps({
            "profile": args.profile, "n_RK": n,
            "maximum_scaled_endpoint_difference": order_difference
        }, sort_keys=True, allow_nan=False), flush=True)
    print("PASS: mixed ANALYTIC/NEWTON group trajectories for N_RK=2/3/4", flush=True)


def order_suite(args):
    final_time = 0.32
    boundary_buffer = 6
    for kind in ("thermal", "drag", "mixed"):
        vertices, q0, initial_t, mode = manufactured(kind)
        # The analytic order comparisons exclude the boundary layer: the
        # outflow-only boundary condition changes the manufactured full-domain
        # solution, while the interior retains the independent uniform-flow
        # reference over this integration interval.
        interior = np.zeros(initial_t.shape, dtype=bool)
        interior[boundary_buffer:-boundary_buffer,
                 boundary_buffer:-boundary_buffer] = True
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
        for n in TEST_STAGES:
            errors = []; fingerprints = []; minima = []; stage_checks = []
            for steps in LEVELS:
                directory = Path(f"order-{kind}-RK{n}-steps{steps}")
                fields, hashes = launch(args.executable, directory, vertices, q0, n, steps,
                                        final_time/steps, drag=int(kind != "thermal"), limiter=0)
                minima.append(admissibility(fields)); fingerprints.append(hashes)
                if np.any(fields["statistics"][:, 0] > fields["statistics"][:, 1]*(1+32768*EPS)):
                    raise AssertionError("order fixture exceeds the retained production CFL bound")
                mass_changes = fields["statistics"][:, 8:10]
                mass_tolerance = 32768*EPS*max(1.0, float(np.sum(fields["q0"][..., 0])))
                if (np.any(mass_changes > mass_tolerance) or
                        np.any(np.diff(mass_changes, axis=0) > mass_tolerance)):
                    raise AssertionError("outflow-only order boundaries increased total mass")
                close("manufactured solid remains zero", fields["q"][..., 4], 0.0)
                close("no raw-stage clipping in order fixture", fields["statistics"][:, 5], 0.0)
                if kind != "thermal":
                    stage_checks.append(implicit_stage_checks(fields, n, final_time/steps))
                    expected_momentum = q0[..., 1:3]/(1+final_time)
                    error = np.max(np.abs(
                        fields["q"][..., 1:3][interior]-
                        expected_momentum[interior]))/RHO
                    if kind == "mixed":
                        error = max(error, float(np.max(np.abs(
                            fields["qp"][..., 3][interior]-expected_t[interior]))))
                else:
                    error = np.max(np.abs(
                        fields["qp"][..., 3][interior]-expected_t[interior]))
                errors.append(float(error))
            rates = np.log2(np.array(errors[:-1])/errors[1:])
            target = 1 if n == 2 else (3 if n == 4 else 2)
            if not np.all((rates[-2:] > target-0.25) & (rates[-2:] < target+0.25)):
                raise AssertionError(f"temporal order {kind} RK{n}: errors={errors}, rates={rates}")
            record = {"profile": args.profile, "kind": kind, "n_RK": n, "expected_order": target,
                      "final_time": final_time, "step_counts": LEVELS,
                      "interior_boundary_buffer_cells": boundary_buffer,
                      "errors": errors,
                      "observed_orders": rates.tolist(), "fingerprints": fingerprints, "admissibility": minima,
                      "implicit_stage_checks": stage_checks,
                      "oracle": "exact constant-coefficient thermal MOL exponential" if kind == "thermal" else
                      "analytic quadratic drag" if kind == "drag" else "independent thermal MOL RK4 plus analytic drag"}
            print("IMEX_ORDER_EVIDENCE "+json.dumps(record, sort_keys=True, allow_nan=False), flush=True)


def stage_count_suite(executable, profile):
    """Require nonzero, diagnostic-specific rejection through both production entry points."""
    generator = Path(__file__).resolve().parent.parent / "TEST_PCCU_LAKE_AT_REST/generate_case.py"
    for mode in ("input", "workspace"):
        for stages in (-1, 0, 1, 2, 5, 6, None):
            if stages is None and mode == "workspace":
                continue  # Omission is a namelist default, not a workspace request.
            directory = Path(f"stage-count-{mode}-{stages}")
            directory.mkdir()
            if mode == "input":
                subprocess.run([sys.executable, str(generator)], cwd=directory, check=True)
                inp = directory / "IMEX_SfloW2D.inp"
                content = inp.read_text()
                assert content.count("N_RK=2,") == 1
                inp.write_text(content.replace("N_RK=2,", "" if stages is None else f"N_RK={stages},"))
            expected = 2 if stages is None else stages
            accepted = expected in (2, 3, 4)
            evidence = []
            for threads in (1, 4):
                env = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE"}
                run = subprocess.run([str(executable), mode, str(expected)], cwd=directory, env=env,
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=30)
                (directory / f"threads-{threads}.log").write_text(run.stdout)
                if accepted:
                    if run.returncode != 0 or "PASS: stage count accepted:" not in run.stdout:
                        raise AssertionError(f"{directory}: valid stage count rejected\n{run.stdout}")
                elif run.returncode == 0 or "FATAL ERROR: N_RK must be 2, 3 or 4" not in run.stdout:
                    raise AssertionError(f"{directory}: missing explicit stage-count rejection\n{run.stdout}")
                evidence.append({"requested_threads": threads, "returncode": run.returncode,
                                 "log_sha256": hashlib.sha256(run.stdout.encode()).hexdigest()})
            print("IMEX_STAGE_COUNT_EVIDENCE "+json.dumps({"profile": profile, "entry": mode,
                  "requested_n_RK": stages, "effective_n_RK": expected, "accepted": accepted,
                  "runs": evidence}, sort_keys=True), flush=True)
    print("PASS: unsupported N_RK rejected by input and workspace; default remains 2", flush=True)
    directory = Path("source-groups-native")
    directory.mkdir()
    subprocess.run([sys.executable, str(generator)], cwd=directory, check=True)
    inp = directory / "IMEX_SfloW2D.inp"
    content = inp.read_text()
    content = content.replace(
        "RADIAL_SOURCE_FLAG=F, COLLAPSING_VOLUME_FLAG=F,",
        "RADIAL_SOURCE_FLAG=F, BOTTOM_RADIAL_SOURCE_FLAG=T, "
        "COLLAPSING_VOLUME_FLAG=F, STOCHASTIC_FLAG=T,",
    )
    content = content.replace("T_START=0.0D0,", "T_START=4.9D0,")
    content = content.replace("T_END=1.0D-1,", "T_END=5.1D0,")
    content = content.replace("DT0=1.0D-4,", "DT0=2.0D-1,")
    inp.write_text(content + """
&STOCHASTIC_PARAMETERS
 std_max=1.0D0,
 tau_stochastic=1.0D0,
 noise_pow_val=1.0D0,
 stochastic_seed=1729,
 stoch_transport_flag=F,
/
&RADIAL_SOURCE_PARAMETERS
 x_source=6.0D0,
 y_source=4.0D0,
 r_source=1.5D0,
 vel_source=1.0D-2,
 T_source=300.0D0,
 time_param=10.0D0,5.0D0,0.0D0,0.0D0,
/
""")
    with (directory / "IMEX_SfloW2D.inp").open("a") as stream:
        stream.write("""
&SOURCE_GROUPS
 source_group_count=3,
 source_group_id(1)=1,
 source_group_order(1)=1,
 source_group_cadence(1)='STAGE',
 source_group_method(1)='NEWTON',
 source_group_term_count(1)=0,
 source_group_target_count(1)=1,
 source_group_targets(1,1)='VELOCITY',
 source_group_id(2)=1,
 source_group_order(2)=1,
 source_group_cadence(2)='STEP',
 source_group_method(2)='EXPLICIT',
 source_group_term_count(2)=1,
 source_group_terms(1,2)='mass_exchange_carrier_liquid_loss',
 source_group_target_count(2)=0,
 source_group_id(3)=2,
 source_group_order(3)=2,
 source_group_cadence(3)='STAGE',
 source_group_method(3)='EXPLICIT',
 source_group_term_count(3)=1,
 source_group_terms(1,3)='bottom_circular_injection',
 source_group_target_count(3)=0
/
""")
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"}
    run = subprocess.run([str(executable), "native-source", "2"], cwd=directory, env=env,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=30)
    (directory / "native-source.log").write_text(run.stdout)
    if run.returncode != 0 or "PASS: stage count accepted:" not in run.stdout or \
            "PASS: production timestep phase order" not in run.stdout:
        raise AssertionError(f"{directory}: native source input flow failed\n{run.stdout}")
    print("PASS: SOURCE_GROUPS read_param flow, production timestep phase order, "
          "and active source-window mass budget", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path); parser.add_argument("profile")
    parser.add_argument("--stage-count-executable", type=Path)
    parser.add_argument("--section", choices=("stages", "order", "analytic", "groups"))
    args = parser.parse_args(); args.executable = args.executable.resolve()
    check_order_conditions()
    if args.stage_count_executable:
        stage_count_suite(args.stage_count_executable.resolve(), args.profile)
    if args.section not in ("order", "analytic", "groups"):
        cases = ("lake_rough", "dynamic_smooth", "dynamic_drainage", "dynamic_compact")
        for name in cases:
            vertices, bed, q = spatial.fixture(name)
            for slope, curvature in ((False, False), (True, True)):
                params = spatial.kt.Parameters(rho_a=1.2, eps_sing=1e-8,
                                               slope_correction=slope, curvature_term=curvature)
                _, diagnostic = spatial.kt.rhs(
                    q, vertices, spatial.DX, spatial.DY, params, return_diagnostics=True)
                dt_x = params.cfl*spatial.DX/float(np.max(np.maximum(
                    -diagnostic["x"]["a_minus"], diagnostic["x"]["a_plus"])))
                dt_y = params.cfl*spatial.DY/float(np.max(np.maximum(
                    -diagnostic["y"]["a_minus"], diagnostic["y"]["a_plus"])))
                bound = min(dt_x, dt_y)
                harmonic = 1.0/(1.0/dt_x+1.0/dt_y)
                for n in TEST_STAGES:
                    for policy, requested, dt in (("retained", 0.0, bound), ("harmonic", harmonic, harmonic)):
                        expected, final, stage_cfl = explicit_stages(q, vertices, params, n, dt)
                        tag = f"{name}-G{int(slope)}-C{int(curvature)}-RK{n}-{policy}"
                        fields, hashes = launch(args.executable, Path(tag), vertices, lift(q), n, 1, requested,
                                                slope=slope, curvature=curvature, limiter=params.limiter)
                        errors = validate_stages(fields, expected, final, dt, bound, name.startswith("lake"))
                        record = {"profile": args.profile, "case": tag, "n_RK": n, "policy": policy,
                                  "dt_used": dt, "dt_directional_minimum": bound, "dt_harmonic": harmonic,
                                  "active_cells": int(fields["active"].sum()), "maximum_scaled_errors": errors,
                                  "stage_CFL_diagnostics": stage_cfl,
                                  "admissibility": admissibility(fields), **hashes}
                        record["intended_final_component_scaled_errors"] = {
                            label: float(np.max(np.abs(fields["raw_final"][..., i]-final[..., i]))) /
                            max(1.0, float(np.max(np.abs(final[..., i]))))
                            for i, label in enumerate(("mass", "mx", "my", "thermal", "solid"))}
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
    if args.section in (None, "order"):
        order_suite(args)
    if args.section in (None, "analytic"):
        analytic_tableau_suite(args)
    if args.section in (None, "groups"):
        mixed_group_tableau_suite(args)


if __name__ == "__main__":
    main()
