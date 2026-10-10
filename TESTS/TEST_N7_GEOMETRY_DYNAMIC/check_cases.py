"""Finite-time 2D geometry gates with raw-stage and full-field IMEX oracles.

The spatial kernels remain the pinned independent Python implementation.
Temporal assembly, final machine-epsilon cleanup and explicitly documented
current-policy cache adapters are independently specified; no production
trajectory is used as a reference field. The historical clean SSPRK2 and
current-policy IMEX comparisons are separate, not interchangeable oracles.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from reference_adapter import (production_dry_momentum_policy, production_slope_geometry,
                               final_dry_face_momenta, stationary_roundoff_policy,
                               auxiliary_tangential_faces, correct_auxiliary_tangential_transport,
                               production_curvature_rhs)
spec = importlib.util.spec_from_file_location(
    "equilibrium_checker", HERE.parent / "TEST_N7_EQUILIBRIUM/check_cases.py")
equilibrium = importlib.util.module_from_spec(spec)
spec.loader.exec_module(equilibrium)
stages = equilibrium.stages
spatial = stages.spatial
core = spatial.core
CONTRACT_PATH = HERE / "contract_v6.json"
CONTRACT = json.loads(CONTRACT_PATH.read_text())
ROUND = CONTRACT["roundoff_epsilon_multiplier"] * np.finfo(float).eps


def compact(state):
    """Map canonical liquid conservative fields to the clean four-field core."""
    return state[..., (0, 1, 2, 4)].copy()


def lift(state):
    """Attach constant-temperature thermal mass without changing the core fields."""
    return stages.lift(state)


def reference(initial, vertices, dx, dy, dt_values, n, slope=False,
              curvature=False, rho_a=1.2, eps_sing=1e-8, collect_bounds=False,
              time_policy=None, stage_observer=None):
    """Advance common initial fields at identical times with independent IMEX weights.

    Reject raw negativity BEFORE the reference kernel's primitive sanitization.
    Record the last known stages and every final dry cleanup separately.
    Optionally return initial directional CFL bounds for each accepted step.
    """
    params = core.SolverParams(rho_a=rho_a, eps_sing=eps_sing,
                               slope_correction=slope, curvature_term=curvature,
                               rest_flux_tol=0.0)
    bed = core.continuous_q1_bed_geometry(vertices)
    sg = production_slope_geometry(bed, dx, dy, core) if slope or curvature else None
    ae, _, be, _ = stages.tableaux(n)
    state = compact(initial)
    minimum = float(state[..., 0].min() / params.rho_c)
    repairs, bounds = [], []
    adapted_cells = 0
    time = 0.0
    dt_old = dt_old_old = time_policy["dt0"] if time_policy else 0.0
    dt_sequence = []
    count = len(dt_values) if time_policy is None else time_policy["maximum_steps"]
    for iteration in range(count):
        if time_policy and time >= time_policy["end"]:
            break
        dt = dt_values[iteration] if time_policy is None else None
        terms, known = [], []
        for i in range(n):
            current = state.copy()
            for weight, term in zip(ae[i, :i], terms):
                current += dt * weight * term
            minimum = min(minimum, float(current[..., 0].min() / params.rho_c))
            if not np.all(np.isfinite(current)) or minimum < -ROUND:
                raise AssertionError("inadmissible raw independent reference stage")
            known.append(lift(current))
            # Retain EVERY raw stage, including stages whose explicit residual
            # is unused in this tableau. No need to evaluate an unused RHS:
            # this is reference-side algebra, not a production optimization.
            if be[i] == 0 and not np.any(ae[i+1:, i]):
                terms.append(np.zeros_like(current))
                continue
            # Production retains positive cell mass below the *face* dry
            # tolerance. Already validated raw cells must not be projected
            # by the prototype's optional global cell sanitizer; face dry
            # handling inside the pinned reconstruction remains unchanged.
            cache = core.prepare_stage_cache2d(current, bed, time, params,
                q_is_sanitized=True, dx=dx, dy=dy, slope_geometry=sg)
            adapted_cells += production_dry_momentum_policy(cache, params, core)
            final_dry_face_momenta(cache, params, core)
            if (collect_bounds or time_policy) and i == 0:
                harmonic, directions = core.max_dt2d(current, bed, dx, dy, time,
                                                     params, stage_cache=cache)
                bounds.append({**directions, "harmonic": harmonic})
                if time_policy:
                    dt = min(directions["dt_x"], directions["dt_y"],
                        time_policy["max_dt"], time_policy["end"]-time,
                        1.1*0.5*(dt_old+dt_old_old))
                    dt_old_old, dt_old = dt_old, dt
            rhs = production_curvature_rhs(current, bed, dx, dy, time, params, cache, core)
            correct_auxiliary_tangential_transport(rhs, cache, dx, dy, params, core)
            stationary_roundoff_policy(rhs, cache, dx, dy, params, core)
            if stage_observer is not None:
                stage_observer(iteration, i, current, cache, rhs, dt)
            terms.append(rhs)
        for weight, term in zip(be, terms):
            state += dt * weight * term
        raw = lift(state)
        minimum = min(minimum, float(state[..., 0].min() / params.rho_c))
        if not np.all(np.isfinite(state)) or minimum < -ROUND:
            raise AssertionError("inadmissible raw independent reference final")
        dry = state[..., 0] < np.finfo(float).eps
        state[dry] = 0
        # User-authorized accepted-state projection; keep raw stages/final
        # assembly observable and retain all positive mass and components.
        unresolved = core.primitive2d(state, params)["h"] <= params.dry_h_tol
        state[..., 1:3][unresolved] = 0
        repairs.append(float(np.max(np.abs(raw-lift(state)))))
        time += dt
        dt_sequence.append(dt)
    if time_policy and abs(time-time_policy["end"]) > ROUND:
        raise AssertionError("independent adaptive trajectory did not reach comparison time")
    return {"q": lift(state), "known": np.array(known), "raw": raw,
            "minimum_raw_h": minimum, "repairs": np.array(repairs), "bounds": bounds,
            "dt_sequence": dt_sequence, "near_dry_momentum_policy_cells": adapted_cells}


def rotate(state, turns=1):
    """Rotate Cartesian fields and momentum consistently on an increasing-y grid."""
    result = state.copy()
    for _ in range(turns):
        result = np.rot90(result, axes=(0, 1)).copy()
        mx = result[..., 1].copy()
        result[..., 1] = result[..., 2]
        result[..., 2] = -mx
    return result


def circular_anisotropy(state, length):
    """Measure the dimensionless fourfold thickness moment, not just axis symmetry.

    Weight with r^4 to regularize the origin: numerator is the Cartesian
    fourth harmonic x^4-6*x^2*y^2+y^4; denominator is h*r^4.
    This detects a D4-symmetric but square-shaped front.
    """
    return abs(circular_fourfold_moment(state, length))


def circular_fourfold_moment(state, length):
    """Return the signed moment so initial sampling and evolution cannot cancel unseen."""
    ny, nx = state.shape[:2]
    x, y = np.meshgrid((np.arange(nx)+0.5)*length/nx-length/2,
                       (np.arange(ny)+0.5)*length/ny-length/2)
    h = state[..., 0] / 1000
    denominator = float(np.sum(h*(x*x+y*y)**2))
    if denominator <= 0:
        raise AssertionError("empty circular isotropy fixture")
    return float(np.sum(h*(x**4-6*x*x*y*y+y**4))) / denominator


def seeded(name, inclined=False):
    """Reuse the original PCG64 beds/states, adding a constant quiet collar.

    The optional locally inclined bed retains exact constant boundary data;
    its compact support exercises a rough continuous-Q1 near-dry front.
    """
    vertices, _, q = spatial.fixture(name)
    if inclined:
        x, y = np.meshgrid(np.arange(spatial.NX+1), np.arange(spatial.NY+1))
        envelope = np.clip((np.minimum(x, spatial.NX-x)-4)/2, 0, 1)
        envelope *= np.clip((np.minimum(y, spatial.NY-y)-4)/2, 0, 1)
        vertices = vertices + envelope*(-0.15*x*spatial.DX+0.08*y*spatial.DY)
    pad = CONTRACT["seeded"]["pad_cells"]
    vertices = np.pad(vertices, pad, mode="edge")
    initial = lift(np.pad(q, ((pad, pad), (pad, pad), (0, 0)), mode="edge"))
    return vertices, initial, spatial.DX, spatial.DY


def mixed(rest=False):
    """Smooth compact quadratic cross term with nonzero u*v*Bxy in wet cells."""
    s = CONTRACT["mixed"]
    x, y = np.meshgrid(np.arange(s["nx"]+1)*s["dx"]-6,
                       np.arange(s["ny"]+1)*s["dy"]-6)
    collar = np.clip((5-np.abs(x))/2, 0, 1)*np.clip((5-np.abs(y))/2, 0, 1)
    vertices = 1 + collar*(0.02*x*y+0.015*x*x+0.01*y*y)
    bed = core.continuous_q1_bed_geometry(vertices)
    xc, yc = np.meshgrid((np.arange(s["nx"])+0.5)*s["dx"]-6,
                         (np.arange(s["ny"])+0.5)*s["dy"]-6)
    support = np.clip((4-np.abs(xc))/2, 0, 1)*np.clip((4-np.abs(yc))/2, 0, 1)
    initial = equilibrium.liquid_state(3-bed.B_center)
    if not rest:
        initial[..., 1] = initial[..., 0]*0.5*support
        initial[..., 2] = initial[..., 0]*0.4*support
    pad = s["pad_cells"]
    return np.pad(vertices, pad, mode="edge"), np.pad(initial,
        ((pad, pad), (pad, pad), (0, 0)), mode="edge"), s["dx"], s["dy"]


def audit_faces(executable, directory, vertices, state, dx, dy, slope, curvature):
    """Observe real sparse reconstruction and compare every active oriented face.

    Use the existing spatial observer without changing its byte format or
    production modules. Check trace positivity/closure and characteristic
    speeds separately; large thermal units cannot hide a momentum defect.
    """
    saved_shape = spatial.NX, spatial.NY
    spatial.NY, spatial.NX = state.shape[:2]
    directory.mkdir(parents=True)
    try:
        spatial.write_fixture(directory/"fixture.inp", vertices, compact(state), slope, curvature, 3)
        # The shared writer uses its spatial fixture spacings; only this
        # generated input record needs the actual spacing of the current case.
        text = (directory/"fixture.inp").read_text().splitlines()
        text[1] = f"{dx:.17g} {dy:.17g}"
        (directory/"fixture.inp").write_text("\n".join(text)+"\n")
        for threads in (1, 4):
            work = directory/f"threads-{threads}"; work.mkdir()
            (work/"fixture.inp").write_bytes((directory/"fixture.inp").read_bytes())
            env = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE"}
            with (work/"solver.log").open("w") as stream:
                run = subprocess.run([str(executable), str(threads)], cwd=work, env=env,
                    stdout=stream, stderr=subprocess.STDOUT, timeout=180)
            if run.returncode:
                raise AssertionError(f"spatial observer failed: {work}\n"+(work/"solver.log").read_text())
        payload = (directory/"threads-1/result.bin").read_bytes()
        if payload != (directory/"threads-4/result.bin").read_bytes():
            raise AssertionError("actual 1/4-thread spatial face payloads differ")
        fields = spatial.read_payload(directory/"threads-1/result.bin")
    finally:
        spatial.NX, spatial.NY = saved_shape
    params = core.SolverParams(rho_a=1.2, eps_sing=1e-8, slope_correction=slope,
                              curvature_term=curvature, rest_flux_tol=0)
    bed = core.continuous_q1_bed_geometry(vertices)
    cache = core.prepare_stage_cache2d(compact(state), bed, 0, params,
        q_is_sanitized=True, dx=dx, dy=dy,
        slope_geometry=production_slope_geometry(bed, dx, dy, core) if slope or curvature else None)
    count = production_dry_momentum_policy(cache, params, core)
    final_dry_face_momenta(cache, params, core)
    active = fields["active"] > 0
    ny, nx = active.shape
    fx = np.zeros((ny, nx+1), dtype=bool); fy = np.zeros((ny+1, nx), dtype=bool)
    fx[:, :-1] |= active; fx[:, 1:] |= active
    fy[:-1] |= active; fy[1:] |= active
    errors = {}
    for normal, mask, entries, sides in (("x", fx, cache["x"], ("L", "R")),
                                        ("y", fy, cache["y"], ("B", "T"))):
        axis = 0 if normal == "x" else 1
        for side, reference_side in zip(sides, ("L", "R")):
            traces = np.stack([entry["cache1d"]["data"][f"q{reference_side}_face"] for entry in entries], axis=axis)
            observed = fields[f"q_{normal}{side}"]
            observed_physical = fields[f"qp_{normal}{side}"]
            if np.min(observed[mask, 0]) < -ROUND or np.min(observed_physical[mask, 0]) < -ROUND:
                raise AssertionError("negative raw reconstructed face")
            normal_component = 1 if normal == "x" else 2
            for pi, ri in ((0, 0), (normal_component, 1), (4, 2)):
                errors[f"{normal}{side}-q{pi}"] = stages.close(
                    "independent near-dry face", observed[..., pi], traces[..., ri], mask=mask)
            tangent_component = 2 if normal == "x" else 1
            tangent_states = compact(state) if normal == "x" else np.swapaxes(compact(state), 0, 1)
            tangent = []
            for row in tangent_states:
                left, right = auxiliary_tangential_faces(row, normal, params, core)
                tangent.append(left if reference_side == "L" else right)
            tangent = np.stack(tangent, axis=axis)
            errors[f"{normal}{side}-tangential"] = stages.close(
                "independent tangential momentum", observed[..., tangent_component],
                traces[..., 0]*tangent, mask=mask)
            eta = np.stack([entry["cache1d"]["data"][f"eta{reference_side}_face"]
                            for entry in entries], axis=axis)
            errors[f"{normal}{side}-eta"] = stages.close(
                "independent free-surface trace", fields[f"eta_{normal}{side}"], eta, mask=mask)
            errors[f"{normal}{side}-thermal"] = stages.close(
                "constant-temperature face closure", observed[..., 3]/spatial.THERMAL_FACTOR,
                traces[..., 0], mask=mask)
            errors[f"{normal}{side}-density"] = stages.close(
                "liquid face density", observed[..., 0]/1000, observed_physical[..., 0], mask=mask)
        for sign, key in (("minus", "a_minus"), ("plus", "a_plus")):
            expected = np.stack([entry["cache1d"][key] for entry in entries], axis=axis)
            errors[f"wave-{normal}-{sign}"] = stages.close(
                "face characteristic bound", fields[f"a{normal}_{sign}"], expected, mask=mask)
    return {"actual_threads": [1, 4], "output_sha256": hashlib.sha256(payload).hexdigest(),
            "maximum_scaled_errors": errors, "active_cells": int(active.sum()),
            "reference_near_dry_policy_cells": count}


def ramp(angle):
    """One-cell continuous ramp at the specified physical angle, with a dry collar."""
    s = CONTRACT["ramp"]; nx = s["nx"]; dx = s["length"]/nx
    x, y = np.meshgrid(np.arange(nx+1)*dx-s["length"]/2,
                       np.arange(nx+1)*dx-s["length"]/2)
    theta = np.deg2rad(angle)
    vertices = s["height"]*np.clip((x*np.cos(theta)+y*np.sin(theta))/dx+0.5, 0, 1)
    bed = core.continuous_q1_bed_geometry(vertices)
    xc, yc = np.meshgrid((np.arange(nx)+0.5)*dx-s["length"]/2,
                         (np.arange(nx)+0.5)*dx-s["length"]/2)
    h = np.maximum(0, 1-(xc*xc+yc*yc)/2.5**2)*(1-bed.B_center)
    initial = equilibrium.liquid_state(h)
    initial[..., 1] = initial[..., 0]*0.2
    initial[..., 2] = initial[..., 0]*0.1
    return vertices, initial, dx, dx


def circular(nx):
    """Smooth radial compact thickness expanding into a genuinely dry exterior."""
    s = CONTRACT["circular"]; dx = s["length"]/nx
    x, y = np.meshgrid((np.arange(nx)+0.5)*dx-s["length"]/2,
                       (np.arange(nx)+0.5)*dx-s["length"]/2)
    h = 0.5*np.maximum(0, 1-(x*x+y*y)/s["radius"]**2)**2
    return np.zeros((nx+1, nx+1)), equilibrium.liquid_state(h), dx, dx


def run(args):
    """Run a fixed complete case inventory and preserve reproducible field evidence."""
    records, differentials, isotropy, failures = [], [], [], []

    def launch(label, fixture, n, count, dt, slope=False, curvature=False, resting=False):
        vertices, initial, dx, dy = fixture
        expected = reference(initial, vertices, dx, dy, [dt]*count, n, slope, curvature)
        fields, fingerprints = stages.launch(args.executable, Path(label), vertices,
            initial, n, count, dt, slope=slope, curvature=curvature, dx=dx, dy=dy)
        diagnostics = equilibrium.validate(fields, initial, count, dry_cleanup=True)
        errors = equilibrium.compare_trajectory(fields, initial, expected["q"],
            expected["known"], count, expected["raw"], expected["repairs"])
        # h and eta use ALL cells, not an active-cell subset or aggregate hash.
        bed = core.continuous_q1_bed_geometry(vertices)
        h = fields["q"][..., 0]/1000
        href = expected["q"][..., 0]/1000
        field = {"h_Linf": float(np.max(np.abs(h-href))),
                 "eta_Linf": float(np.max(np.abs((h+bed.B_center)-(href+bed.B_center))))}
        if max(field.values()) > ROUND*count*max(1, float(initial[..., 0].max()/1000)):
            raise AssertionError("full thickness/free-surface comparison failed")
        record = {"case": label, "n_RK": n, "steps": count, "dt": dt,
                  "time": dt*count, "slope": slope, "curvature": curvature,
                  "shape": list(initial.shape[:2]), "fingerprints": fingerprints,
                  "admissibility": diagnostics, "reference_errors": errors,
                  "full_fields": field, "reference_minimum_raw_h": expected["minimum_raw_h"]}
        record["reference_near_dry_momentum_policy_cells"] = expected["near_dry_momentum_policy_cells"]
        if args.spatial_executable:
            record["last_step_stage_face_audits"] = [audit_faces(args.spatial_executable,
                Path(label)/f"face-stage-{i+1}", vertices, fields["solved"][i], dx, dy, slope, curvature)
                for i in range(n)]
        if resting:
            record["equilibrium_errors"] = equilibrium.field_error(
                "mixed lake unchanged", fields["q"], initial, initial, ROUND)
        records.append(record)
        print("PASS:", label, flush=True)
        return fields["q"], record

    s = CONTRACT["seeded"]
    for name in s["fixtures"] + ["inclined_compact"]:
        data = seeded("dynamic_compact" if name == "inclined_compact" else name,
                      name == "inclined_compact")
        for slope, curvature in s["flags"]:
            for n in CONTRACT["stages"]:
                launch(f"{name}-S{int(slope)}C{int(curvature)}-RK{n}", data,
                       n, s["steps"], s["dt"], slope, curvature)

    s = CONTRACT["mixed"]
    data = mixed()
    sg = core.build_slope_geometry_2d(core.continuous_q1_bed_geometry(data[0]), data[2], data[3])
    # Fixture guard prevents the curvature differential becoming a vacuous toggle.
    u = data[1][..., 1]/data[1][..., 0]; v = data[1][..., 2]/data[1][..., 0]
    if np.max(np.abs(2*sg.Bxy_center*u*v)) < 1e-4:
        raise AssertionError("mixed Hessian fixture has no active cross contribution")
    for n in CONTRACT["stages"]:
        results = {}
        for slope, curvature in s["flags"]:
            for rest in (False, True):
                q, record = launch(f"mixed-{'lake' if rest else 'dynamic'}-S{int(slope)}C{int(curvature)}-RK{n}",
                    mixed(rest), n, s["steps"], s["dt"], slope, curvature, rest)
                if not rest:
                    results[slope, curvature] = q
        for slope in (False, True):
            delta = float(np.max(np.abs(results[slope, True]-results[slope, False])/
                                 equilibrium.scales(data[1])))
            if delta < s["minimum_curvature_differential"]:
                raise AssertionError("dynamic curvature toggle has no measurable effect")
            differentials.append({"n_RK": n, "slope": slope, "maximum_scaled_differential": delta})

    s = CONTRACT["ramp"]
    for angle in s["angles_degrees"]:
        data = ramp(angle)
        for n in CONTRACT["stages"]:
            first = None
            for turn in s["rotations"]:
                vertices, initial, dx, dy = data
                q, record = launch(f"ramp-angle{angle}-rotation{turn}-RK{n}",
                    (np.rot90(vertices, turn).copy(), rotate(initial, turn), dx, dy),
                    n, s["steps"], s["dt"])
                if turn == 0:
                    first = q
                else:
                    record["rotation_errors"] = equilibrium.field_error(
                        "quarter-turn covariance", q, rotate(first, turn), initial, ROUND*s["steps"])

    s = CONTRACT["circular"]; count = round(s["time"]/s["dt"])
    for n in CONTRACT["stages"]:
        metrics, evolution_metrics = [], []
        for nx in s["nx"]:
            data = circular(nx)
            q, record = launch(f"circular-nx{nx}-RK{n}", data, n, count, s["dt"])
            record["rotation_errors"] = equilibrium.field_error(
                "circular quarter-turn symmetry", q, rotate(q), data[1], ROUND*count)
            record["anisotropy"] = circular_anisotropy(q, s["length"])
            record["initial_anisotropy"] = circular_anisotropy(data[1], s["length"])
            record["signed_initial_fourfold_moment"] = circular_fourfold_moment(data[1], s["length"])
            record["signed_final_fourfold_moment"] = circular_fourfold_moment(q, s["length"])
            record["evolution_fourfold_moment_change"] = (
                record["signed_final_fourfold_moment"]-record["signed_initial_fourfold_moment"])
            if record["anisotropy"] > s["maximum_fourfold_anisotropy"]:
                raise AssertionError("circular grid-aligned distortion exceeds frozen limit")
            metrics.append(record["anisotropy"])
            evolution_metrics.append(abs(record["evolution_fourfold_moment_change"]))
        ratio = metrics[-1]/max(metrics[0], np.finfo(float).tiny)
        evolution_ratio = evolution_metrics[-1]/max(evolution_metrics[0], np.finfo(float).tiny)
        if evolution_ratio > s["maximum_finest_over_coarsest_evolution_anisotropy"]:
            failures.append(f"RK{n}: circular evolution distortion does not decrease with refinement ({evolution_ratio:.9g})")
        isotropy.append({"n_RK": n, "anisotropy": metrics,
                        "historical_final_only_finest_over_coarsest": ratio,
                        "evolution_anisotropy": evolution_metrics,
                        "evolution_finest_over_coarsest": evolution_ratio})

    if len(records) != 69:
        raise AssertionError("incomplete N7-C fixed inventory")
    result = {"profile": args.profile, "cases": records, "case_count": len(records),
        "actual_threads": [1, 4], "solver_runs": 2*len(records),
        "contract_sha256": hashlib.sha256(CONTRACT_PATH.read_bytes()).hexdigest(),
        "curvature_differentials": differentials, "isotropy_refinements": isotropy,
        "reference_core_sha256": spatial.CORE_HASHES, "CFL_policy_changed": False,
        "status": "failed" if failures else "passed", "failures": failures}
    Path("evidence.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    if failures:
        raise AssertionError("; ".join(failures))
    print(f"PASS: N7-C {args.profile}, {len(records)} finite-time cases", flush=True)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("executable", type=Path)
    cli.add_argument("profile", choices=CONTRACT["profiles"])
    cli.add_argument("--spatial-executable", type=Path)
    run(cli.parse_args())
