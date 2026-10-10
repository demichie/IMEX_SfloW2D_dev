"""Frozen N7-B multicomponent equilibrium, transport and raw budget gates.

The scalar contact oracle never calls production or the one-solid prototype.
Thermodynamic identities are independently expressed as specific-volume and
heat-capacity sums; no clipped primitive field is a positivity oracle.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "TEST_IMEX_STAGES"))
import compare_stages as stages

CONTRACT = json.loads((HERE / "contract.json").read_text())
P = CONTRACT["properties"]
ROUND = CONTRACT["roundoff_epsilon_multiplier"] * np.finfo(float).eps
RHO_A = P["pressure"] / (P["R_a"] * P["temperature"])


def composition(name, profile, neutral=True):
    """Build nonconstant fractions; neutral changes lie in the specific-volume nullspace."""
    count = {"liquid": 2, "gas": 3, "gas-liquid": 4}[name]
    base = np.array([0.15, 0.10] if name == "liquid" else [0.20, 0.15, 0.05] +
                    ([0.20] if name == "gas-liquid" else []))
    delta = np.zeros(count)
    delta[0] = 0.025
    if neutral:
        if name == "liquid":
            delta[1] = -delta[0]*(1/P["rho_s"][0]-1/P["rho_l"])/(1/P["rho_s"][1]-1/P["rho_l"])
        else:
            air_volume = P["R_a"]*P["temperature"]/P["pressure"]
            delta[1] = -0.015
            if name == "gas-liquid":
                delta[3] = 0.02
            numerator = np.dot(delta[:2], 1/np.array(P["rho_s"])-air_volume)
            if name == "gas-liquid":
                numerator += delta[3]*(1/P["rho_l"]-air_volume)
            delta[2] = -numerator / ((P["R_g"]-P["R_a"])*P["temperature"]/P["pressure"])
    return base + profile[..., None]*delta, base, delta


def properties(name, fractions, temperature):
    """Independent current liquid/solid, gas/solid and gas/liquid/solid closures."""
    remainder = 1-fractions.sum(axis=-1)
    cp = fractions[..., :2] @ np.array(P["cp_s"])
    volume = fractions[..., :2] @ (1/np.array(P["rho_s"]))
    if name == "liquid":
        cp += remainder*P["cp_l"]
        volume += remainder/P["rho_l"]
    else:
        cp += remainder*P["cp_a"] + fractions[..., 2]*P["cp_g"]
        volume += (remainder*P["R_a"]+fractions[..., 2]*P["R_g"])*temperature/P["pressure"]
        if name == "gas-liquid":
            cp += fractions[..., 3]*P["cp_l"]
            volume += fractions[..., 3]/P["rho_l"]
    return 1/volume, cp, remainder


def conservative(name, fractions, h, velocity=0.0, axis="x"):
    """Pack physical cell averages into the canonical thermal conservative state."""
    rho, cp, _ = properties(name, fractions, P["temperature"])
    mass = rho*h
    result = np.zeros(fractions.shape[:-1]+(4+fractions.shape[-1],))
    result[..., 0] = mass
    result[..., 1 if axis == "x" else 2] = mass*velocity
    result[..., 3] = mass*cp*P["temperature"]
    result[..., 4:] = mass[..., None]*fractions
    return result


def decode(name, state):
    """Decode raw states without clipping mass, fractions, temperature or velocity."""
    if not np.all(np.isfinite(state)) or np.any(state[..., 0] <= 0):
        raise AssertionError("nonfinite or nonpositive raw wet state")
    fractions = state[..., 4:]/state[..., :1]
    _, cp, remainder = properties(name, fractions, P["temperature"])
    temperature = state[..., 3]/(state[..., 0]*cp)
    rho, _, _ = properties(name, fractions, temperature)
    return fractions, temperature, rho, state[..., 0]/rho, remainder


def field_scales(initial):
    """Scale each component independently; thermal units never mask momentum errors."""
    scale = np.maximum(1, np.max(np.abs(initial), axis=(0, 1)))
    scale[1:3] = scale[0]*np.sqrt(P["gravity"]*max(1, float(np.max(initial[..., 0]))/P["rho_l"]))
    return scale


def compare(label, observed, expected, initial, limit):
    """Reject any nonfinite or excessive component-wise full-field discrepancy."""
    if observed.shape != expected.shape or not np.all(np.isfinite(observed)) or not np.all(np.isfinite(expected)):
        raise AssertionError(label+": invalid field")
    error = np.max(np.abs(observed-expected)/field_scales(initial), axis=tuple(range(observed.ndim-1)))
    if np.any(error > limit):
        raise AssertionError(f"{label}: component errors {error}, frozen limit {limit}")
    return error.tolist()


def pulse_cell_means(nx, shift=0.0):
    """Independent analytic translated cos^4 pulse, averaged with 32-point quadrature."""
    settings = CONTRACT["contact"]
    nodes, weights = np.polynomial.legendre.leggauss(32)
    dx = settings["length"]/nx
    x = (np.arange(nx)+0.5)*dx
    distance = (x[:, None]+0.5*dx*nodes-shift-settings["profile_center"])/settings["profile_half_width"]
    value = np.where(np.abs(distance) < 1, np.cos(0.5*np.pi*distance)**4, 0)
    return 0.5*(value@weights)


def scalar_rhs(value, dx, velocity, sound):
    """Independent generalized-minmod scalar CU transport on a uniform state."""
    slope = np.zeros_like(value)
    a, b = value[2:]-value[1:-1], value[1:-1]-value[:-2]
    c = 0.5*(value[2:]-value[:-2])
    theta = CONTRACT["contact"]["theta"]
    slope[1:-1] = np.where((a*b > 0) & (a*c > 0),
                           np.sign(c)*np.minimum(np.abs(c), theta*np.minimum(np.abs(a), np.abs(b))), 0)
    minus, plus = value-0.5*slope, value+0.5*slope
    left, right = np.r_[minus[0], plus], np.r_[minus, plus[-1]]
    a_plus, a_minus = max(0, velocity+sound), min(0, velocity-sound)
    flux = (a_plus*velocity*left-a_minus*velocity*right+a_plus*a_minus*(right-left))/(a_plus-a_minus)
    return -np.diff(flux)/dx


def scalar_reference(profile, dx, velocity, sound, dt, count, n):
    """Apply independently specified current IMEX explicit weights, without a prototype stepper."""
    ae, _, be, _ = stages.tableaux(n)
    value = profile.copy()
    for _ in range(count):
        terms, known = [], []
        for i in range(n):
            stage = value.copy()
            for weight, rhs in zip(ae[i, :i], terms):
                stage += dt*weight*rhs
            known.append(stage)
            terms.append(scalar_rhs(stage, dx, velocity, sound))
        for weight, rhs in zip(be, terms):
            value += dt*weight*rhs
    return value, np.array(known)


def read_payload(path, nx, ny, nv, n, steps):
    """Read exactly the unchanged base observer format with the selected variable extent."""
    raw = np.fromfile(path, dtype=np.float64)
    cursor, fields = 0, {}
    for name, shape in (("q0", (ny, nx, nv)), ("q", (ny, nx, nv)), ("qp", (ny, nx, nv+2)),
                        ("active", (ny, nx)), ("known", (n, ny, nx, nv)), ("solved", (n, ny, nx, nv)),
                        ("raw_final", (ny, nx, nv)), ("implicit_statuses", (n, ny, nx)), ("statistics", (steps, 15))):
        count = int(np.prod(shape))
        fields[name] = raw[cursor:cursor+count].reshape(shape)
        cursor += count
    if cursor != raw.size or not np.all(np.isfinite(raw)):
        raise AssertionError("invalid or nonfinite base payload")
    raw = np.fromfile(path.with_name("composition.bin"), dtype=np.float64)
    cursor = 0
    for name, shape in (("composition", (steps, 2*(nv+1)+5)),
                        ("W", (ny, nx, nv+2)), ("E", (ny, nx, nv+2)),
                        ("S", (ny, nx, nv+2)), ("N", (ny, nx, nv+2)),
                        ("L", (ny, nx+1, nv)), ("R", (ny, nx+1, nv)),
                        ("B", (ny+1, nx, nv)), ("T", (ny+1, nx, nv)), ("rhs", (ny, nx, nv))):
        count = int(np.prod(shape))
        fields[name] = raw[cursor:cursor+count].reshape(shape)
        cursor += count
    if cursor != raw.size or not np.all(np.isfinite(raw)):
        raise AssertionError("invalid or nonfinite composition payload")
    return fields


def validate(name, fields, initial, count):
    """Check all-step raw budgets/face extrema and independently decode retained raw fields."""
    statistics, audits = fields["statistics"], fields["composition"]
    nv = initial.shape[-1]
    if statistics.shape != (count, 15) or audits.shape != (count, 2*(nv+1)+5):
        raise AssertionError("incomplete all-step composition diagnostics")
    inventory = initial.sum(axis=(0, 1))
    carrier = inventory[0]-inventory[4:].sum()
    scales = np.maximum(1, np.abs(np.r_[inventory, carrier]))
    conservative_indices = [0, 3, *range(4, nv), nv]
    budget = np.maximum(audits[:, :nv+1], audits[:, nv+1:2*(nv+1)])/scales
    if np.any(budget[:, conservative_indices] > ROUND):
        raise AssertionError("raw/final component or residual-carrier budget")
    if np.any(statistics[:, 2:5] < -ROUND) or np.any(statistics[:, 6] != 0):
        raise AssertionError("negative raw mass or failed local solve")
    if np.any(statistics[:, 12] < 273-ROUND) or np.any(statistics[:, 13:15] < -ROUND):
        raise AssertionError("all-step raw temperature/component/carrier inadmissible")
    if np.any(statistics[:, 0] > statistics[:, 1]*(1+ROUND)):
        raise AssertionError("diagnostic step exceeds retained CFL")
    if np.any(statistics[:, 5] > ROUND*field_scales(initial).max()):
        raise AssertionError("unexpected final repair in fully wet fixture")
    extra = audits[:, 2*(nv+1):]
    if np.any(extra[:, [0, 3, 4]] > ROUND):
        raise AssertionError("raw reconstructed face positivity/TVD/temperature")
    if np.any(extra[:, 1] > ROUND*max(1, float(np.max(initial[..., 0])))):
        raise AssertionError("closed boundary leaks mass")
    minima = {"component": 1.0, "carrier": 1.0, "temperature": float("inf"), "h": float("inf")}
    for key in ("q0", "q", "known", "solved", "raw_final", "L", "R", "B", "T"):
        fractions, temperature, rho, h, remainder = decode(name, fields[key])
        if np.min(fractions) < -ROUND or np.min(remainder) < -ROUND or np.min(temperature) < 273-ROUND:
            raise AssertionError(key+": independent raw admissibility")
        minima["component"] = min(minima["component"], float(fractions.min()))
        minima["carrier"] = min(minima["carrier"], float(remainder.min()))
        minima["temperature"] = min(minima["temperature"], float(temperature.min()))
        minima["h"] = min(minima["h"], float(h.min()))
    # Independently reconstruct final conservative endpoints from PRE-closure
    # cell-owned traces, testing the thermodynamic closure in every component.
    face_errors = {}
    for primitive, conservative_key, slicer, axis in (
            ("W", "R", (slice(None), slice(None, -1)), "x"),
            ("E", "L", (slice(None), slice(1, None)), "x"),
            ("S", "T", (slice(None, -1), slice(None)), "y"),
            ("N", "B", (slice(1, None), slice(None)), "y")):
        trace = fields[primitive]
        rho, cp, remainder = properties(name, trace[..., 4:nv], trace[..., 3])
        expected = np.zeros(trace.shape[:-1]+(nv,))
        expected[..., 0] = rho*trace[..., 0]
        expected[..., 1:3] = rho[..., None]*trace[..., 1:3]
        expected[..., 3] = expected[..., 0]*cp*trace[..., 3]
        expected[..., 4:] = expected[..., :1]*trace[..., 4:nv]
        face_errors[primitive] = compare("independent face closure "+primitive, fields[conservative_key][slicer],
                                         expected, initial, ROUND)
    return {"max_component_budget": np.max(budget[:, conservative_indices], axis=0).tolist(),
            "budget_components": conservative_indices, "raw_minima": minima,
            "max_preclosure_face_violation": float(extra[:, 0].max()),
            "max_boundary_mass_flux": float(extra[:, 1].max()),
            "max_momentum_residual": float(extra[:, 2].max()),
            "max_dt_over_CFL": float(np.max(statistics[:, 0]/statistics[:, 1])),
            "face_closure_errors": face_errors}


def launch(executable, label, name, initial, vertices, n, steps, dt, dx, dy, flags, limiter):
    """Require actual one/four teams and bitwise equality of BOTH raw payloads."""
    folder = Path(label); folder.mkdir()
    stages.write_fixture(folder/"fixture.inp", vertices, initial, *flags, limiter, dt, dx, dy)
    mode = {"liquid": 1, "gas": 2, "gas-liquid": 3}[name]
    outputs = []
    for team in (1, 4):
        work = folder/f"threads-{team}"; work.mkdir()
        (work/"fixture.inp").write_bytes((folder/"fixture.inp").read_bytes())
        env = {**os.environ, "OMP_NUM_THREADS": str(team), "OMP_DYNAMIC": "FALSE"}
        with (work/"solver.log").open("w") as log:
            run = subprocess.run([str(executable), str(team), str(n), str(steps), "0", "1", str(mode)],
                                 cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180)
        if run.returncode:
            raise AssertionError(f"{work}: exit {run.returncode}\n"+(work/"solver.log").read_text())
        ny, nx, nv = initial.shape
        outputs.append(read_payload(work/"result.bin", nx, ny, nv, n, steps))
    hashes = {}
    for file in ("result.bin", "composition.bin"):
        a, b = (folder/f"threads-{team}"/file for team in (1, 4))
        if a.read_bytes() != b.read_bytes():
            raise AssertionError(label+": actual 1/4-thread payload mismatch: "+file)
        hashes[file] = hashlib.sha256(a.read_bytes()).hexdigest()
    hashes["input"] = hashlib.sha256((folder/"fixture.inp").read_bytes()).hexdigest()
    return outputs[0], {"sha256": hashes, "actual_threads": [1, 4]}


def run(args):
    """Execute only frozen cases, independently checking equilibrium and signed transport."""
    records, refinements = [], []
    eq = CONTRACT["equilibrium"]
    for name in CONTRACT["closures"]:
        for dimension in (1, 2):
            nx, ny = eq["nx"], eq["ny_2d"] if dimension == 2 else 1
            x, y = np.meshgrid((np.arange(nx)+0.5)*eq["spacing"], (np.arange(ny)+0.5)*eq["spacing"])
            profile = 0.5+0.2*np.sin(0.6*x)+0.2*np.cos(0.7*y)
            for family in eq["families"]:
                fractions, base, delta = composition(name, profile, neutral=family=="constant-density-Q1")
                rho, _, _ = properties(name, fractions, P["temperature"])
                if family == "constant-density-Q1":
                    if np.ptp(rho) > ROUND*float(rho.max()):
                        raise AssertionError("independent density-nullspace construction failed")
                    xv, yv = np.meshgrid(np.arange(nx+1)*eq["spacing"], np.arange(ny+1)*eq["spacing"])
                    vertices = 1+0.1*np.cos(0.4*xv)
                    if dimension == 2:
                        vertices += 0.2*np.sin(0.5*yv)
                    bed = 0.25*(vertices[:-1, :-1]+vertices[1:, :-1]+vertices[:-1, 1:]+vertices[1:, 1:])
                    h = eq["eta"]-bed
                    flags_list, limiter = eq["flags"], eq["density_limiter"]
                else:
                    vertices = np.zeros((ny+1, nx+1))
                    rho0, _, _ = properties(name, base, P["temperature"])
                    h = np.sqrt((rho0-RHO_A)/(rho-RHO_A))
                    if np.ptp(h) < 1e-5:
                        raise AssertionError("pressure equilibrium must have nonconstant eta and Gamma")
                    pressure = P["gravity"]*(rho-RHO_A)*h*h
                    if np.ptp(pressure) > ROUND*float(pressure.max()):
                        raise AssertionError("independent hydrostatic pressure identity failed")
                    flags_list, limiter = eq["pressure_flat_flags"], eq["pressure_flat_limiter"]
                initial = conservative(name, fractions, h)
                for flags in flags_list:
                    for n in CONTRACT["stages"]:
                        label = f"lake-{family}-{name}-{dimension}D-S{int(flags[0])}C{int(flags[1])}-RK{n}"
                        fields, fingerprints = launch(args.executable, label, name, initial, vertices, n, eq["steps"],
                                                      eq["dt"], eq["spacing"], eq["spacing"], flags, limiter)
                        diagnostics = validate(name, fields, initial, eq["steps"])
                        errors = {key: compare(label+": equilibrium "+key, fields[key],
                                              np.broadcast_to(initial, fields[key].shape), initial, ROUND)
                                  for key in ("q", "known", "solved", "raw_final")}
                        force_scale = max(1, float(np.max(initial[..., 0])))*P["gravity"]/eq["spacing"]
                        residual = diagnostics["max_momentum_residual"]/force_scale
                        if residual > ROUND:
                            raise AssertionError("equilibrium residual exceeds frozen force-scaled limit")
                        records.append({"case": label, "steps": eq["steps"], "fingerprints": fingerprints,
                                        "diagnostics": diagnostics, "equilibrium_errors": errors,
                                        "force_scaled_residual": residual})
                        print("PASS:", label, flush=True)
    if args.equilibrium_only:
        if len(records) != 90:
            raise AssertionError("incomplete frozen equilibrium inventory")
        Path("equilibrium_evidence.json").write_text(json.dumps({
            "status": "equilibrium subset passed; transport and N7-B closure remain open",
            "profile": args.profile, "contract": CONTRACT, "cases": records,
        }, indent=2, allow_nan=False)+"\n")
        print(f"PASS: {args.profile} N7-B equilibrium subset, 90 cases, actual 1/4 threads", flush=True)
        return
    settings = CONTRACT["contact"]
    for name in CONTRACT["closures"]:
        for axis in settings["axes"]:
            for sign in settings["directions"]:
                for n in CONTRACT["stages"]:
                    series = []
                    for count in (settings["nx"] if axis == "x" else settings["y_refinements"]):
                        dx = settings["length"]/count
                        profile = pulse_cell_means(count)
                        fractions, base, delta = composition(name, profile)
                        rho0, _, _ = properties(name, base, P["temperature"])
                        sound = np.sqrt(P["gravity"]*(1-RHO_A/rho0))
                        coordinates = (np.arange(count)+0.5)*dx
                        rise = np.clip((coordinates-5)/5, 0, 1)
                        fall = np.clip((35-coordinates)/5, 0, 1)
                        velocity = sign*settings["velocity"]*(rise*rise*(3-2*rise))*(fall*fall*(3-2*fall))
                        fractions = fractions[None] if axis == "x" else fractions[:, None]
                        vel = velocity[None] if axis == "x" else velocity[:, None]
                        initial = conservative(name, fractions, np.ones(vel.shape), vel, axis)
                        nx, ny = initial.shape[1], initial.shape[0]
                        vertices = np.zeros((ny+1, nx+1))
                        steps = int(np.ceil(settings["time"]/(settings["dt_dx_factor"]*dx)))
                        dt = settings["time"]/steps
                        label = f"contact-{name}-{axis}-sign{sign}-nx{count}-RK{n}"
                        fields, fingerprints = launch(args.executable, label, name, initial, vertices, n, steps,
                                                      dt, dx, dx, (False, False), settings["limiter"])
                        diagnostics = validate(name, fields, initial, steps)
                        expected_profile, known_profile = scalar_reference(profile, dx, sign*settings["velocity"],
                                                                           sound, dt, steps, n)
                        interior = (coordinates >= settings["comparison_interval"][0]) & (coordinates <= settings["comparison_interval"][1])
                        mask = interior[None] if axis == "x" else interior[:, None]

                        def lift(values):
                            values = values[None] if axis == "x" else values[:, None]
                            ys, _, _ = composition(name, values)
                            return conservative(name, ys, np.ones(values.shape), sign*settings["velocity"], axis)

                        expected = lift(expected_profile)
                        known = np.array([lift(value) for value in known_profile])
                        errors = {key: compare("independent scalar contact "+key,
                                               fields[key][:, mask] if fields[key].ndim == 4 else fields[key][mask],
                                               known[:, mask] if fields[key].ndim == 4 else expected[mask],
                                               initial, ROUND*steps)
                                  for key in ("q", "known", "solved", "raw_final")}
                        ys, _, _, _, _ = decode(name, fields["q"])
                        observed_profile = (ys[..., 0]-base[0])/delta[0]
                        observed_profile = observed_profile[0] if axis == "x" else observed_profile[:, 0]
                        analytic = pulse_cell_means(count, sign*settings["velocity"]*settings["time"])
                        l1 = float(np.mean(np.abs(observed_profile[interior]-analytic[interior])))
                        linf = float(np.max(np.abs(observed_profile[interior]-analytic[interior])))
                        if l1 > settings["maximum_profile_L1"] or linf > settings["maximum_profile_Linf"]:
                            raise AssertionError("analytic contact error exceeds frozen limit")
                        mean0 = float(np.sum(coordinates*profile)/np.sum(profile))
                        mean1 = float(np.sum(coordinates*observed_profile)/np.sum(observed_profile))
                        # A nontrivial pulse must move in the prescribed sign;
                        # the independent scalar field gate controls its magnitude.
                        if sign*(mean1-mean0) <= 0:
                            raise AssertionError("material contact did not move in the prescribed direction")
                        series.append(l1)
                        records.append({"case": label, "steps": steps, "fingerprints": fingerprints,
                                        "diagnostics": diagnostics, "reference_errors": errors,
                                        "analytic_L1": l1, "analytic_Linf": linf, "centroid_displacement": mean1-mean0})
                        print("PASS:", label, flush=True)
                    if axis == "x":
                        ratios = (np.array(series[1:])/series[:-1]).tolist()
                        if max(ratios) >= settings["maximum_refinement_ratio"]:
                            raise AssertionError("contact does not satisfy frozen refinement criterion")
                        refinements.append({"closure": name, "sign": sign, "n_RK": n, "L1": series, "ratios": ratios})
    expected_count = 90+72
    if len(records) != expected_count:
        raise AssertionError("incomplete frozen N7-B inventory")
    evidence = {"profile": args.profile, "contract": CONTRACT, "cases": records, "refinements": refinements}
    Path("evidence.json").write_text(json.dumps(evidence, indent=2, allow_nan=False)+"\n")
    print(f"PASS: {args.profile} N7-B, {len(records)} cases, actual 1/4 threads", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("profile", choices=CONTRACT["profiles"])
    parser.add_argument("--equilibrium-only", action="store_true",
                        help="Run all 90 frozen equilibria without claiming transport or N7-B closure")
    args = parser.parse_args()
    try:
        run(args)
    except (AssertionError, subprocess.SubprocessError) as error:
        # A failed gate must remain a failure. Persist its frozen contract and
        # cause before exiting; this is never a partial acceptance certificate.
        Path("failure.json").write_text(json.dumps({
            "status": "failed; N7-B remains open",
            "profile": args.profile,
            "error": str(error),
            "expected_inventory": 162,
            "contract_sha256": hashlib.sha256((HERE/"contract.json").read_bytes()).hexdigest(),
            "contract": CONTRACT,
        }, indent=2, allow_nan=False)+"\n")
        raise
