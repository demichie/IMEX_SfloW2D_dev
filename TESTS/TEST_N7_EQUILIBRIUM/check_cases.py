"""N7-A finite-time production gates with criteria frozen in contract.json.

The existing driver observes raw known/solved/final states on every step and
checks the actual OpenMP team. This harness never uses primitive clipping as
a positivity oracle. Reference kernels are checksum-pinned by compare_reference.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

TEST_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TEST_DIR.parent / "TEST_IMEX_STAGES"))
import compare_stages as stages  # noqa: E402
import hp_pccu_1d_core as core1  # noqa: E402

CONTRACT = json.loads((TEST_DIR / "contract.json").read_text())
ROUND = CONTRACT["roundoff_epsilon_multiplier"] * np.finfo(float).eps
RHO, CP, TEMP, GRAV, RHO_A = 1000.0, 4180.0, 300.0, 9.81, 1.2


def scales(initial):
    """Keep thermal-energy units from masking a mass or momentum defect."""
    mass = max(1.0, float(np.max(np.abs(initial[..., 0]))))
    return np.array([mass, mass * np.sqrt(GRAV * max(1.0, mass / RHO)),
                     mass * np.sqrt(GRAV * max(1.0, mass / RHO)),
                     max(1.0, float(np.max(np.abs(initial[..., 3])))), mass])


def field_error(label, actual, expected, initial, limit):
    """Compare each conservative component at every cell against a fixed bound."""
    if not np.all(np.isfinite(actual)) or not np.all(np.isfinite(expected)):
        raise AssertionError(f"{label}: nonfinite field")
    error = np.max(np.abs(actual - expected) / scales(initial), axis=tuple(range(actual.ndim - 1)))
    if np.any(error > limit):
        raise AssertionError(f"{label}: component errors {error}, limit {limit}")
    return error.tolist()


def validate(fields, initial, steps, budget=True, dry_cleanup=False):
    """Reject raw negativity, failed solves, hidden repairs and CFL violations."""
    diagnostics = stages.admissibility(fields)
    statistics = fields["statistics"]
    if np.any(statistics[:, 0] > statistics[:, 1] * (1 + ROUND)):
        raise AssertionError("prescribed diagnostic step exceeds production CFL")
    if not dry_cleanup and diagnostics["maximum_final_repair"] > ROUND * scales(initial).max():
        raise AssertionError("final repair hides an inadmissible raw state")
    total = max(1.0, float(np.sum(initial[..., 0])))
    if budget and np.max(np.abs(statistics[:, 8:10])) > ROUND * total:
        raise AssertionError("raw/final mass budget drift")
    if abs(statistics[-1, 11] - np.sum(statistics[:, 0])) > ROUND * steps:
        raise AssertionError("wrong finite integration time")
    diagnostics["maximum_relative_mass_drift"] = float(np.max(np.abs(statistics[:, 8:10])) / total)
    diagnostics["maximum_dt_over_CFL"] = float(np.max(statistics[:, 0] / statistics[:, 1]))
    return diagnostics


def liquid_state(h, u=0.0, temperature=TEMP):
    """Create canonical thermal conservative states for a pure liquid fixture."""
    state = np.zeros(np.shape(h) + (5,))
    state[..., 0] = RHO * h
    state[..., 1] = RHO * h * u
    state[..., 3] = RHO * h * CP * temperature
    return state


def lake_fixture(name, dimension):
    """Sample authoritative shared Q1 vertices; fill using their exact cell mean."""
    settings = CONTRACT["lake"]
    nx, ny = settings["nx"], settings["ny_2d"] if dimension == 2 else 1
    dx, dy = settings["dx"], settings["dy"]
    x, y = np.meshgrid(np.arange(nx + 1) * dx, np.arange(ny + 1) * dy)
    if dimension == 1:
        y = np.zeros_like(y)
    if name == "smooth":
        vertices = 1 + 0.15 * np.cos(0.7 * x) + 0.1 * np.sin(0.6 * y)
    elif name == "one-cell":
        vertices = 1 + 1.2 * np.clip((x - 4) / dx, 0, 1)
        if dimension == 2:
            vertices += 0.8 * np.clip((y - 3) / dy, 0, 1)
    elif name == "parabolic":
        vertices = 1 + 0.03 * (x - 4)**2
        if dimension == 2:
            vertices += 0.02 * (y - 3)**2
    else:
        raise AssertionError("unlisted lake fixture")
    bed = stages.spatial.core.continuous_q1_bed_geometry(vertices)
    state = liquid_state(settings["eta"] - bed.B_center)
    if state[..., 0].min() <= 0 or vertices.max() >= settings["eta"]:
        raise AssertionError("lake must be fully wet including all Q1 faces")
    return vertices, state, dx, dy


def compact(state):
    """Map production zero-transverse states into the frozen 1D core layout."""
    return state[0, :, (0, 1, 4)].T.copy()


def lift(state):
    """Reinsert transverse momentum and the passively transported thermal mass."""
    result = np.zeros((1, state.shape[0], 5))
    result[0, :, 0:2] = state[:, :2]
    result[0, :, 4] = state[:, 2]
    result[0, :, 3] = CP * TEMP * state[:, 0]
    return result


def reference_1d(initial, vertices, dx, dt, count, n, slope=False, curvature=False):
    """Advance unchanged spatial kernels using independently specified IMEX weights."""
    params = core1.SolverParams(rho_a=RHO_A, eps_sing=1e-8,
                               slope_correction=slope, curvature_term=curvature)
    bed = core1.continuous_bed_geometry(vertices[0])
    ae, _, be, _ = stages.tableaux(n)
    state = compact(initial)
    minimum = float(core1.raw_thickness(state, params).min())
    dry_repairs = []
    for iteration in range(count):
        rhs = []
        known = []
        for i in range(n):
            stage = state.copy()
            for weight, term in zip(ae[i, :i], rhs):
                stage += dt * weight * term
            minimum = min(minimum, float(core1.raw_thickness(stage, params).min()))
            if minimum < -ROUND:
                raise AssertionError("negative raw reference stage before kernel sanitization")
            known.append(lift(stage))
            rhs.append(core1.rhs_spatial(stage, bed, dx, iteration * dt, params,
                                        core1.transmissive_boundary))
        for weight, term in zip(be, rhs):
            state += dt * weight * term
        minimum = min(minimum, float(core1.raw_thickness(state, params).min()))
        if minimum < -ROUND:
            raise AssertionError("negative raw reference final before next kernel sanitization")
        raw = lift(state)
        # Production's existing final canonicalization discards momentum only
        # when total mass is below machine epsilon (not the interface h cutoff).
        # Capture and compare raw states BEFORE this explicitly audited mapping.
        dry = state[:, 0] < np.finfo(float).eps
        dry_repairs.append(float(np.max(np.abs(raw[0, dry]))) if np.any(dry) else 0.0)
        state[dry] = 0
    return lift(state), np.array(known), minimum, raw, np.array(dry_repairs)


def thermal_reference(initial, dx, dt, count, n, velocity):
    """Independent first-order central-upwind scalar transport on a uniform liquid."""
    temperature = initial[0, :, 3] / (RHO * CP)
    speed = np.sqrt(GRAV * (1 - RHO_A / RHO))
    ae, _, be, _ = stages.tableaux(n)

    def rhs(value):
        left = np.r_[value[0], value]
        right = np.r_[value, value[-1]]
        flux = 0.5 * (velocity + speed) * left + 0.5 * (velocity - speed) * right
        return -np.diff(flux) / dx

    for _ in range(count):
        terms, known = [], []
        for i in range(n):
            value = temperature.copy()
            for weight, term in zip(ae[i, :i], terms):
                value += dt * weight * term
            known.append(liquid_state(np.ones((1, value.size)), velocity, value))
            terms.append(rhs(value))
        for weight, term in zip(be, terms):
            temperature += dt * weight * term
    return liquid_state(np.ones((1, temperature.size)), velocity, temperature), np.array(known)


def compare_trajectory(fields, initial, final, known, count, raw=None, dry_repairs=None):
    """Check full final arrays and all raw states of the last retained step."""
    limit = ROUND * count
    errors = {"final": field_error("finite-time final", fields["q"], final, initial, limit)}
    active = fields["active"] > 0
    # The last active workset can grow at a dry front; compare final everywhere,
    # but never read the observer's intentionally zero inactive entries as states.
    for key in ("known", "solved"):
        observed = fields[key]
        observed_active = np.any(observed[..., 0] != 0, axis=0) | active
        errors[key] = field_error(key, observed[:, observed_active], known[:, observed_active], initial, limit)
    observed_active = np.any(fields["known"][..., 0] != 0, axis=0) | active
    errors["raw_final"] = field_error("raw_final", fields["raw_final"][observed_active],
                                      (final if raw is None else raw)[observed_active], initial, limit)
    if dry_repairs is not None:
        # A dry-only momentum reset is allowed, but any unpredicted repair fails.
        errors["dry_cleanup"] = float(np.max(np.abs(fields["statistics"][:, 5] - dry_repairs)) / scales(initial)[:3].max())
        if errors["dry_cleanup"] > limit:
            raise AssertionError("production repair differs from independently predicted dry cleanup")
    return errors


def ritter_cell_means(nx, length, time):
    """Integrate the analytic rarefaction polynomial exactly over each cell."""
    c = np.sqrt(GRAV * (1 - RHO_A / RHO))
    edges = np.linspace(0, length, nx + 1) - CONTRACT["ritter"]["dam_x"]
    left, right = edges[:-1], edges[1:]
    a, b = np.maximum(left, -c * time), np.minimum(right, 2 * c * time)
    b = np.maximum(a, b)
    wet_left = np.maximum(0, np.minimum(right, -c * time) - left)

    def fh(x):
        return (4*c*c*x - 2*c*x*x/time + x**3/(3*time*time)) / (9*c*c)

    def fm(x):
        return (4*c**3*x + 0*x*x - c*x**3/time**2 + x**4/(4*time**3)) * 2 / (27*c*c)

    return (wet_left + fh(b) - fh(a)) / (length / nx), (fm(b) - fm(a)) / (length / nx)


def run(args):
    """Execute exactly the frozen inventory and emit machine-readable evidence."""
    records, refinements = [], []

    def launch(label, vertices, initial, n, count, dt, dx, dy=1.0,
               slope=False, curvature=False, limiter=3, budget=True, dry_cleanup=False):
        fields, hashes = stages.launch(args.executable, Path(label), vertices, initial,
                                       n, count, dt, slope=slope, curvature=curvature,
                                       limiter=limiter, dx=dx, dy=dy)
        record = {"case": label, "n_RK": n, "steps": count, "dt": dt,
                  "shape": list(initial.shape[:2]), "slope": slope, "curvature": curvature,
                  "fingerprints": hashes, "admissibility": validate(fields, initial, count, budget, dry_cleanup)}
        records.append(record)
        return fields, record

    for bed_name in CONTRACT["lake"]["beds"]:
        for dimension in (1, 2):
            vertices, initial, dx, dy = lake_fixture(bed_name, dimension)
            for slope, curvature in CONTRACT["lake"]["flags"]:
                params = stages.spatial.core.SolverParams(rho_a=RHO_A, eps_sing=1e-8,
                                                          slope_correction=slope, curvature_term=curvature)
                bed = stages.spatial.core.continuous_q1_bed_geometry(vertices)
                reference_q = np.zeros(initial.shape[:2] + (4,))
                reference_q[..., :3] = initial[..., :3]
                residual = stages.spatial.core.rhs2d(reference_q, bed, dx, dy, 0, params)
                force_scale = RHO * GRAV * (initial[..., 0].max() / RHO)**2 / dx
                if np.max(np.abs(residual)) > ROUND * force_scale:
                    raise AssertionError("independent lake residual is not zero")
                for n in CONTRACT["stages"]:
                    count, dt = CONTRACT["lake"]["steps"], CONTRACT["lake"]["dt"]
                    fields, record = launch(f"lake-{bed_name}-{dimension}D-S{int(slope)}C{int(curvature)}-RK{n}",
                                            vertices, initial, n, count, dt, dx, dy, slope, curvature)
                    record["equilibrium_errors"] = field_error("unchanged lake", fields["q"], initial, initial, ROUND)
                    for key in ("known", "solved", "raw_final"):
                        field_error(key, fields[key], initial, initial, ROUND)
                    # Stage 2's known value is q_start + dt*a21*R(q_start).
                    ae = stages.tableaux(n)[0]
                    inferred = (fields["known"][1] - fields["known"][0]) / (dt * ae[1, 0])
                    record["production_lake_residual_scaled"] = float(np.max(np.abs(inferred[..., :3])) / force_scale)
                    record["reference_lake_residual_scaled"] = float(np.max(np.abs(residual)) / force_scale)
                    if record["production_lake_residual_scaled"] > ROUND:
                        raise AssertionError("nonzero production lake residual")
                    print("PASS:", record["case"], flush=True)

    for kind in ("advection", "ritter"):
        settings = CONTRACT[kind]
        for n in CONTRACT["stages"]:
            errors = []
            for nx in settings["nx"]:
                dx = settings["length"] / nx
                x = (np.arange(nx) + 0.5) * dx
                vertices = np.zeros((2, nx + 1))
                count = int(np.ceil(settings["time"] / (settings["dt_dx_factor"] * dx)))
                dt = settings["time"] / count
                if kind == "advection":
                    temperature = settings["temperature_base"] + settings["temperature_amplitude"] * np.exp(
                        -((x - settings["center"]) / settings["width"])**2)
                    initial = liquid_state(np.ones((1, nx)), settings["velocity"], temperature)
                    expected, known = thermal_reference(initial, dx, dt, count, n, settings["velocity"])
                    limiter = settings["limiter"]
                else:
                    initial = liquid_state((x < settings["dam_x"])[None, :] * settings["left_h"])
                    expected, known, raw_minimum, raw, dry_repairs = reference_1d(initial, vertices, dx, dt, count, n)
                    limiter = 3
                fields, record = launch(f"{kind}-nx{nx}-RK{n}", vertices, initial, n, count, dt, dx,
                                        limiter=limiter, budget=kind != "advection", dry_cleanup=kind == "ritter")
                record["reference_errors"] = compare_trajectory(fields, initial, expected, known, count,
                                                                  raw if kind == "ritter" else None,
                                                                  dry_repairs if kind == "ritter" else None)
                if kind == "advection":
                    exact = settings["temperature_base"] + settings["temperature_amplitude"] * np.exp(
                        -((x - settings["center"] - settings["velocity"] * settings["time"]) / settings["width"])**2)
                    actual = fields["q"][0, :, 3] / (CP * fields["q"][0, :, 0])
                    error = float(np.mean(np.abs(actual - exact)))
                    linf = float(np.max(np.abs(actual - exact)))
                    if error > settings["maximum_L1_temperature_error"] or linf > settings["maximum_Linf_temperature_error"]:
                        raise AssertionError("advection analytic error exceeds frozen limit")
                    # Equal inlet/outlet mass of the uniform moving liquid must cancel.
                    if record["admissibility"]["maximum_relative_mass_drift"] > ROUND:
                        raise AssertionError("uniform advection mass budget drift")
                    record["analytic_Linf"] = linf
                else:
                    h, hu = ritter_cell_means(nx, settings["length"], settings["time"])
                    error = float(np.mean(np.abs(fields["q"][0, :, 0] / RHO - h)))
                    momentum_error = float(np.mean(np.abs(fields["q"][0, :, 1] / RHO - hu)) / np.sqrt(GRAV))
                    if error > settings["maximum_mean_h_error"] or momentum_error > settings["maximum_mean_scaled_momentum_error"]:
                        raise AssertionError("Ritter analytic error exceeds frozen limit")
                    record["analytic_momentum_L1"] = momentum_error
                    record["reference_raw_minimum_h"] = raw_minimum
                record["analytic_L1"] = error
                errors.append(error)
                print("PASS:", record["case"], "L1=", error, flush=True)
            ratios = (np.array(errors[1:]) / errors[:-1]).tolist()
            if not all(0 < ratio < settings["maximum_refinement_ratio"] for ratio in ratios):
                raise AssertionError(f"{kind} RK{n}: nonconverging refinement {ratios}")
            refinements.append({"case": kind, "n_RK": n, "L1": errors, "ratios": ratios})

    settings = CONTRACT["excavation"]
    nx = settings["nx"]; dx = settings["length"] / nx
    x = np.arange(nx + 1) * dx
    west, east = settings["walls"]
    weight = np.clip(np.minimum((x - west) / dx, (east - x) / dx), 0, 1)
    base = 10 - np.tan(np.deg2rad(settings["slope_degrees"])) * x
    vertices = np.tile(base - settings["depth"] * weight, (2, 1))
    initial = liquid_state((settings["depth"] * 0.5 * (weight[:-1] + weight[1:]))[None, :])
    for slope, curvature in CONTRACT["lake"]["flags"]:
        for n in CONTRACT["stages"]:
            count, dt = settings["steps"], settings["dt"]
            expected, known, raw_minimum, raw, dry_repairs = reference_1d(initial, vertices, dx, dt, count, n, slope, curvature)
            fields, record = launch(f"excavation-S{int(slope)}C{int(curvature)}-RK{n}",
                                    vertices, initial, n, count, dt, dx, slope=slope, curvature=curvature, dry_cleanup=True)
            record["reference_errors"] = compare_trajectory(fields, initial, expected, known, count, raw, dry_repairs)
            uphill = (np.arange(nx) + 0.5) * dx < west
            uphill_mass = float(np.sum(np.abs(fields["q"][0, uphill, 0])))
            record["relative_uphill_mass"] = uphill_mass / np.sum(initial[..., 0])
            record["reference_raw_minimum_h"] = raw_minimum
            if record["relative_uphill_mass"] > ROUND:
                raise AssertionError("resolved uphill excavation leakage")
            if fields["q"][..., 1].sum() <= 0:
                raise AssertionError("excavation did not accelerate downslope")
            print("PASS:", record["case"], flush=True)

    if len(records) != 102:
        raise AssertionError("incomplete N7-A fixture inventory")
    result = {"profile": args.profile, "contract_sha256": hashlib.sha256((TEST_DIR / "contract.json").read_bytes()).hexdigest(),
              "cases": records, "refinements": refinements, "case_count": len(records),
              "solver_runs": 2 * len(records), "actual_threads": [1, 4]}
    Path("evidence.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"PASS: N7-A {args.profile}, {len(records)} cases, actual 1/4-thread full payloads identical", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("profile", choices=CONTRACT["profiles"])
    run(parser.parse_args())
