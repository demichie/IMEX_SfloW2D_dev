"""Compare complete production spatial arrays with the pinned HP-PCCU cores.

The explicit Python core and retained IMEX code have different 2-D timestep
combination policies. Directional bounds are strict comparisons; the combined
policy difference is quantified, never relabeled as reference agreement.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

REFERENCE = Path(__file__).parent / "reference"
CORE_HASHES = {
    "hp_pccu_1d_core.py": "ffcb695ff8d8d2bd06e433d70b2a519c0b0012d959eb664b617deea83e8f64e2",
    "hp_pccu_2d_core.py": "c63331366cadf9f04d802a3f058088b9c175399b2d6dc7015068a0ad110f1577",
}
for filename, checksum in CORE_HASHES.items():
    if hashlib.sha256((REFERENCE / filename).read_bytes()).hexdigest() != checksum:
        raise AssertionError(f"unapproved reference core: {filename}")
sys.path.insert(0, str(REFERENCE))
import hp_pccu_2d_core as core  # noqa: E402

NX, NY = 18, 16
DX, DY = 0.75, 1.25
THERMAL_FACTOR = 300.0 * 4180.0
EPS = np.finfo(np.float64).eps


def fixture(name):
    """Deterministic Q1 beds/states with a constant collar for common boundary semantics."""
    rng = np.random.Generator(np.random.PCG64(20261009))
    x, y = np.meshgrid(np.arange(NX + 1), np.arange(NY + 1))
    vx = np.clip((np.minimum(x, NX - x) - 4) / 2.0, 0, 1)
    vy = np.clip((np.minimum(y, NY - y) - 4) / 2.0, 0, 1)
    envelope = vx * vy
    vertices = np.ones((NY + 1, NX + 1))
    if "smooth" in name:
        vertices += envelope * (0.25*np.sin(0.4*x) + 0.15*np.cos(0.5*y) + 0.02*x*y)
    elif "rough" in name or "drainage" in name or "compact" in name:
        vertices += 2.0 * envelope * rng.uniform(-1, 1, vertices.shape)
    elif "ramp" in name:
        vertices += 4.0 * envelope * ((x >= 9) & (y >= 7))
    bed = core.continuous_q1_bed_geometry(vertices)
    xc, yc = np.meshgrid(np.arange(NX) + 0.5, np.arange(NY) + 0.5)
    ec = np.clip((np.minimum(xc, NX-xc)-4)/2.0, 0, 1) * np.clip((np.minimum(yc, NY-yc)-4)/2.0, 0, 1)
    if name.startswith("lake"):
        h = np.max(bed.B_center) + 2.0 - bed.B_center
        u = np.zeros_like(h); v = np.zeros_like(h)
    else:
        h = 0.7 + ec*(0.12*np.sin(0.7*xc) + 0.08*np.cos(0.4*yc))
        u = ec*(0.7*np.sin(0.4*xc) + 0.2*np.cos(0.3*yc))
        v = ec*(-0.5*np.cos(0.6*yc) + 0.15*np.sin(0.2*xc))
        if name == "dynamic_drainage":
            h = 0.7 + ec*(0.003 + 0.002*rng.random(h.shape) - 0.7)
        elif name == "dynamic_compact":
            support = (xc-NX/2)**2 + (yc-NY/2)**2 < 1.9**2
            h = support * (0.3 + 0.15*rng.random(h.shape))
    q = core.state2d_from_h_alpha_uv(h, np.zeros_like(h), u, v,
                                  core.SolverParams(rho_a=1.2, eps_sing=1e-8))
    return vertices, bed, q


def write_fixture(path, vertices, q, slope, curvature, limiter):
    """Write generated numeric input, not a production solver input schema."""
    production = np.zeros((NY, NX, 5))
    production[..., :3] = q[..., :3]
    production[..., 3] = THERMAL_FACTOR*q[..., 0]
    production[..., 4] = q[..., 3]
    with path.open("w") as stream:
        stream.write(f"{NX} {NY} {int(slope)} {int(curvature)} {limiter}\n{DX:.17g} {DY:.17g}\n")
        np.savetxt(stream, vertices.reshape(-1, 1), fmt="%.17g")
        np.savetxt(stream, production.reshape(-1, 1), fmt="%.17g")


def read_payload(path):
    """Decode version-one test stream; Fortran variable-first is a C-last NumPy view."""
    raw = np.fromfile(path, dtype=np.float64)
    fields = {}
    cursor = 0

    def take(name, shape):
        nonlocal cursor
        count = int(np.prod(shape))
        fields[name] = raw[cursor:cursor+count].reshape(shape)
        cursor += count

    take("dt", (1,))
    for name, components in (("q", 5), ("qp", 7), ("spatial", 5), ("sources", 5)):
        take(name, (NY, NX, components))
    take("active", (NY, NX))
    take("vertices", (NY+1, NX+1))
    take("bed", (NY, NX)); take("bed_x", (NY, NX+1)); take("bed_y", (NY+1, NX))
    for name in ("bx", "by", "bxx", "bxy", "byy", "G"):
        take(name, (NY, NX))
    take("Gx", (NY, NX+1)); take("Gy", (NY+1, NX))
    for name, shape in (("ax_minus", (NY, NX+1)), ("ax_plus", (NY, NX+1)),
                        ("ay_minus", (NY+1, NX)), ("ay_plus", (NY+1, NX))):
        take(name, shape)
    for name, shape in (("xL", (NY, NX+1)), ("xR", (NY, NX+1)),
                        ("yB", (NY+1, NX)), ("yT", (NY+1, NX))):
        take("q_"+name, shape+(5,))
    for name, shape in (("xL", (NY, NX+1)), ("xR", (NY, NX+1)),
                        ("yB", (NY+1, NX)), ("yT", (NY+1, NX))):
        take("qp_"+name, shape+(7,))
    for name, shape in (("xL", (NY, NX+1)), ("xR", (NY, NX+1)),
                        ("yB", (NY+1, NX)), ("yT", (NY+1, NX))):
        take("eta_"+name, shape)
    if cursor != len(raw) or not np.all(np.isfinite(raw)):
        raise AssertionError("incorrect stream length or nonfinite payload")
    return fields


def compare(fields, q, bed, params, name):
    """Compare every active cell/face and full geometric fields, not aggregate fingerprints."""
    cache = core.prepare_stage_cache2d(q, bed, 0.0, params, dx=DX, dy=DY)
    expected_rhs, diagnostics = core.rhs2d(q, bed, DX, DY, 0.0, params,
                                         return_diagnostics=True, stage_cache=cache)
    expected_dt, dt_diagnostics = core.max_dt2d(q, bed, DX, DY, 0.0, params, stage_cache=cache)
    active = fields["active"] > 0
    if name == "dynamic_compact" and (np.all(active) or not np.any(active & (fields["qp"][..., 0] == 0))):
        raise AssertionError("compact fixture did not exercise sparse and dry active worksets")
    fx = np.zeros((NY, NX+1), dtype=bool); fy = np.zeros((NY+1, NX), dtype=bool)
    fx[:, :-1] |= active; fx[:, 1:] |= active
    fy[:-1, :] |= active; fy[1:, :] |= active
    errors = {}

    def check(label, actual, expected, mask=None, scale=None):
        actual, expected = np.broadcast_arrays(np.asarray(actual), np.asarray(expected))
        if mask is not None:
            actual = actual[mask]; expected = expected[mask]
        if not actual.size or not np.all(np.isfinite(actual)) or not np.all(np.isfinite(expected)):
            raise AssertionError(f"{name}: invalid comparison {label}")
        magnitude = max(1.0, float(np.max(np.abs(expected)))) if scale is None else max(1.0, scale)
        error = float(np.max(np.abs(actual-expected)))
        errors[label] = error/magnitude
        if error > 8192*EPS*magnitude:
            index = np.unravel_index(int(np.argmax(np.abs(actual-expected))), actual.shape)
            raise AssertionError(f"{name}: {label}, error={error:.17g}, limit={8192*EPS*magnitude:.17g}, "
                                 f"index={index}, actual={actual[index]}, expected={expected[index]}")

    check("bed_centers", fields["bed"], bed.B_center)
    check("bed_x_faces", fields["bed_x"], bed.Bx_face)
    check("bed_y_faces", fields["bed_y"], bed.By_face)
    if params.slope_correction or params.curvature_term:
        sg = core.build_slope_geometry_2d(bed, DX, DY)
        for label, attr in (("bx", "Bx_center"), ("by", "By_center"), ("bxx", "Bxx_center"),
                            ("bxy", "Bxy_center"), ("byy", "Byy_center")):
            check("geometry_"+label, fields[label], getattr(sg, attr))
        expected_g = sg.G_center if params.slope_correction else np.ones((NY, NX))
        expected_gx = sg.Gx_face if params.slope_correction else np.ones((NY, NX+1))
        expected_gy = sg.Gy_face if params.slope_correction else np.ones((NY+1, NX))
    else:
        expected_g = np.ones((NY, NX)); expected_gx = np.ones((NY, NX+1)); expected_gy = np.ones((NY+1, NX))
    check("G_center", fields["G"], expected_g)
    check("G_x_faces", fields["Gx"], expected_gx)
    check("G_y_faces", fields["Gy"], expected_gy)

    # The inactive contract is explicit: no residual is read there in Fortran.
    # Require the independent full-grid reference to be zero outside the active
    # workset in this fixture, so comparison cannot hide an omitted active flux.
    if np.any(~active):
        check("inactive_reference_rhs", expected_rhs[~active], 0.0)
    rhs = -fields["spatial"] + fields["sources"]
    hmax = float(np.max(core.raw_thickness2d(q, params)))
    lake_scale = params.rho_c*params.g*max(1.0, hmax*hmax)/min(DX, DY)
    for production_index, reference_index, label in ((0, 0, "mass"), (1, 1, "momentum_x"),
                                                     (2, 2, "momentum_y"), (4, 3, "solid")):
        check("residual_"+label, rhs[..., production_index], expected_rhs[..., reference_index], active,
              lake_scale if name.startswith("lake") else None)
    check("thermal_transport", fields["spatial"][..., 3]/THERMAL_FACTOR, fields["spatial"][..., 0], active)
    check("no_thermal_source", fields["sources"][..., 3], 0.0, active)
    check("no_mass_source", fields["sources"][..., 0], 0.0, active)
    check("closed_collar_mass_budget", np.sum(rhs[..., 0])*DX*DY, 0.0,
          scale=float(np.sum(np.abs(rhs[..., 0]))*DX*DY))
    check("curvature_x", fields["sources"][..., 1], diagnostics["curvature_source_x"], active)
    check("curvature_y", fields["sources"][..., 2], diagnostics["curvature_source_y"], active)
    if name.startswith("lake"):
        check("lake_residual_zero", rhs[..., :3], 0.0, active, lake_scale)

    for normal, mask, entries, sides in (("x", fx, cache["x"], ("L", "R")),
                                          ("y", fy, cache["y"], ("B", "T"))):
        axis = 0 if normal == "x" else 1
        for label, attr in (("minus", "a_minus"), ("plus", "a_plus")):
            expected = np.stack([entry["cache1d"][attr] for entry in entries], axis=axis)
            check(f"wave_{normal}_{label}", fields[f"a{normal}_{label}"], expected, mask)
        for side, reference_side in zip(sides, ("L", "R")):
            compact = np.stack([entry["cache1d"]["data"][f"q{reference_side}_face"] for entry in entries], axis=axis)
            expected_eta = np.stack([entry["cache1d"]["data"][f"eta{reference_side}_face"] for entry in entries], axis=axis)
            actual_q = fields[f"q_{normal}{side}"]
            actual_qp = fields[f"qp_{normal}{side}"]
            if np.any(actual_qp[..., 0][mask] < 0.0):
                raise AssertionError("negative final production trace thickness")
            normal_component = 1 if normal == "x" else 2
            tangent_component = 2 if normal == "x" else 1
            for pi, ri in ((0, 0), (normal_component, 1), (4, 2)):
                check(f"trace_{normal}{side}_q{pi}", actual_q[..., pi], compact[..., ri], mask)
            tangent_states = q if normal == "x" else np.swapaxes(q, 0, 1)
            tangent = []
            for row in tangent_states:
                vt = core._tangential_velocity_slice(row, normal, params)
                vt_left, vt_right = core._reconstruct_scalar_faces(vt, params)
                tangent.append(vt_left if reference_side == "L" else vt_right)
            tangent = np.stack(tangent, axis=axis)
            check(f"trace_{normal}{side}_tangential", actual_q[..., tangent_component], compact[..., 0]*tangent, mask)
            check(f"trace_{normal}{side}_thermal", actual_q[..., 3]/THERMAL_FACTOR, compact[..., 0], mask)
            check(f"trace_{normal}{side}_eta", fields[f"eta_{normal}{side}"], expected_eta, mask)
            expected_h = (compact[..., 0]-compact[..., 2])/params.rho_c + compact[..., 2]/params.rho_s
            expected_normal = core.desingularized_velocity(compact[..., 0],compact[..., 1],params.eps_sing)
            expected_tangent = core.desingularized_velocity(compact[..., 0],compact[..., 0]*tangent,params.eps_sing)
            check(f"primitive_{normal}{side}_h", actual_qp[..., 0], expected_h, mask)
            check(f"primitive_{normal}{side}_normal_velocity", actual_qp[..., normal_component+4], expected_normal, mask)
            check(f"primitive_{normal}{side}_tangent_velocity", actual_qp[..., tangent_component+4], expected_tangent, mask)
            check(f"primitive_{normal}{side}_temperature", actual_qp[..., 3], 300.0, mask)
            check(f"primitive_{normal}{side}_normal_volume_momentum", actual_qp[..., normal_component], compact[..., 1]/params.rho_c, mask)

    check("normal_x_momentum_mean", 0.5*(fields["q_xL"][:, 1:, 1]+fields["q_xR"][:, :-1, 1]), q[..., 1], active)
    check("normal_y_momentum_mean", 0.5*(fields["q_yB"][1:, :, 2]+fields["q_yT"][:-1, :, 2]), q[..., 2], active)
    check("x_mass_mean", 0.5*(fields["q_xL"][:, 1:, 0]+fields["q_xR"][:, :-1, 0]), q[..., 0], active)
    check("y_mass_mean", 0.5*(fields["q_yB"][1:, :, 0]+fields["q_yT"][:-1, :, 0]), q[..., 0], active)

    dt_x = params.cfl*DX/float(np.max(np.maximum(-fields["ax_minus"], fields["ax_plus"])))
    dt_y = params.cfl*DY/float(np.max(np.maximum(-fields["ay_minus"], fields["ay_plus"])))
    check("directional_dt_x", dt_x, dt_diagnostics["dt_x"])
    check("directional_dt_y", dt_y, dt_diagnostics["dt_y"])
    actual_dt = float(fields["dt"][0])
    check("retained_IMEX_dt_policy", actual_dt, min(dt_x, dt_y))
    policy_ratio = actual_dt/float(expected_dt)
    if not 1.0 < policy_ratio <= 2.0+8192*EPS:
        raise AssertionError("unexpected CFL combination-policy relationship")
    return {
        "maximum_scaled_errors": errors, "active_cells": int(np.sum(active)),
        "dry_active_cells": int(np.sum(active & (fields["qp"][..., 0] == 0))),
        "inactive_cells": int(np.sum(~active)), "dt_IMEX": actual_dt,
        "dt_reference_SSPRK2": float(expected_dt), "dt_x": dt_x, "dt_y": dt_y,
        "timestep_ratio_IMEX_over_reference": policy_ratio,
        "combined_CFL_policy_matches_reference": False,
        "minimum_G": float(np.min(fields["G"])),
        "minimum_raw_Euler_mass_at_IMEX_dt": float(np.min(fields["q"][..., 0]+actual_dt*rhs[..., 0])),
        "scope": "Fixed one-component liquid closure, shared constant boundary collar, no IMEX advancement.",
    }


def negative_controls(fields, q, bed, params):
    """Require diagnostic failures for deliberate residual/trace/geometry/speed/CFL errors."""
    mutations = (
        ("residual_mass", "spatial", (8, 9, 0), 1.0),
        ("trace_xL_tangential", "q_xL", (8, 9, 2), 1.0),
        ("primitive_xL_normal_velocity", "qp_xL", (8, 9, 5), 0.1),
        ("G_x_faces", "Gx", (8, 9), 0.1),
        ("wave_x_plus", "ax_plus", (8, 9), 0.1),
        ("retained_IMEX_dt_policy", "dt", (0,), 0.01),
        ("inactive_reference_rhs", "active", (8, 9), -1.0),
    )
    for label, field, index, amount in mutations:
        bad = {key: value.copy() for key, value in fields.items()}
        bad[field][index] += amount
        try:
            compare(bad, q, bed, params, "dynamic_rough")
        except AssertionError as error:
            if label not in str(error):
                raise AssertionError(f"negative control {label} failed for a different reason: {error}") from error
        else:
            raise AssertionError(f"negative control did not detect {label}")
    print(f"PASS: {len(mutations)} spatial comparator negative controls", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("profile")
    parser.add_argument("--case", choices=("lake_flat", "lake_smooth", "lake_rough", "lake_ramp", "dynamic_smooth",
                                           "dynamic_rough", "dynamic_drainage", "dynamic_compact"))
    parser.add_argument("--require-reference-timestep", action="store_true",
                        help="fail unless the full combined 2-D timestep matches the explicit reference")
    args = parser.parse_args()
    cases = [args.case] if args.case else ["lake_flat", "lake_smooth", "lake_rough", "lake_ramp", "dynamic_smooth",
                                         "dynamic_rough", "dynamic_drainage", "dynamic_compact"]
    payloads = {}
    for name in cases:
        vertices, bed, q = fixture(name)
        for slope, curvature in ((False, False), (True, False), (False, True), (True, True)):
            tag = f"{name}-G{int(slope)}-C{int(curvature)}"
            params = core.SolverParams(rho_a=1.2, eps_sing=1e-8, slope_correction=slope,
                                       curvature_term=curvature)
            results = []
            paths = []
            for threads in (1, 4):
                directory = Path(tag)/f"threads-{threads}"
                directory.mkdir(parents=True)
                write_fixture(directory/"fixture.inp", vertices, q, slope, curvature, params.limiter)
                environment = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE"}
                with (directory/"solver.log").open("w") as log:
                    run = subprocess.run([str(args.executable), str(threads)], cwd=directory,
                                         env=environment, stdout=log, stderr=subprocess.STDOUT, timeout=90)
                if run.returncode:
                    raise AssertionError(f"{tag}, threads={threads}: exit {run.returncode}\n"
                                         +(directory/"solver.log").read_text())
                paths.append(directory/"result.bin")
                record = compare(read_payload(paths[-1]), q, bed, params, name)
                record.update(case=tag, profile=args.profile, actual_threads=threads,
                              reference_core_sha256=CORE_HASHES,
                              input_sha256=hashlib.sha256((directory/"fixture.inp").read_bytes()).hexdigest(),
                              output_sha256=hashlib.sha256(paths[-1].read_bytes()).hexdigest())
                results.append(record)
                if threads == 1:
                    payloads[(name, slope, curvature)] = read_payload(paths[-1])
            if paths[0].read_bytes() != paths[1].read_bytes():
                raise AssertionError(f"{tag}: production 1/4-thread payloads differ")
            print(f"PASS: {args.profile} {tag}, complete residual/geometry/traces/directional CFL, actual 1/4 threads", flush=True)
            for record in results:
                print("N7_SPATIAL_EVIDENCE "+json.dumps(record, sort_keys=True, allow_nan=False), flush=True)
            if name == "dynamic_rough" and slope and not curvature:
                negative_controls(payloads[(name, slope, curvature)], q, bed, params)
            if args.require_reference_timestep:
                raise AssertionError(f"combined timestep differs: IMEX={results[0]['dt_IMEX']}, "
                                     f"reference SSPRK2={results[0]['dt_reference_SSPRK2']}")
    # Curvature is a separate source: toggling it cannot alter transport,
    # reconstruction, geometry or characteristic bounds at fixed slope flag.
    for name in cases:
        for slope in (False, True):
            before = payloads[(name, slope, False)]
            after = payloads[(name, slope, True)]
            for key in before:
                if key == "sources":
                    continue
                if before[key].tobytes() != after[key].tobytes():
                    raise AssertionError(f"curvature toggle changed spatial field {name}, G={slope}, {key}")
            # Curvature arithmetic can return -0 instead of initialized +0;
            # both are exact zero sources, not a physical equilibrium change.
            if name.startswith("lake") and not np.array_equal(before["sources"],after["sources"]):
                raise AssertionError("curvature changed lake-at-rest sources")
    print("PASS: curvature is separate; toggling it leaves spatial fields and CFL bitwise unchanged", flush=True)
    print("CFL_POLICY_DIFFERENCE: retained IMEX min(dt_x,dt_y) versus explicit-reference harmonic combination", flush=True)


if __name__ == "__main__":
    main()
