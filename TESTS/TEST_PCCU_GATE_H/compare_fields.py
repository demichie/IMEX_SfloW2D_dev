"""Reproducible all-cell Gate-H comparison, distinct from historical port norms."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "geometry_dynamic_checker", HERE.parent / "TEST_N7_GEOMETRY_DYNAMIC/check_cases.py")
dynamic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dynamic)
CONTRACT_PATH = HERE/"field_contract_v3.json"
CONTRACT = json.loads(CONTRACT_PATH.read_text())
if CONTRACT["reference_core_sha256"] != dynamic.spatial.CORE_HASHES:
    raise AssertionError("field contract does not pin the approved clean kernels")


def compare(actual, expected, bed, initial_h, steps):
    """Reject any-cell h/eta disagreement beyond the pre-frozen roundoff bound."""
    if actual.shape != expected.shape or bed.shape != actual.shape:
        raise AssertionError("full-field grid shape mismatch")
    if not all(np.all(np.isfinite(f)) for f in (actual, expected, bed)):
        raise AssertionError("nonfinite full-field comparison")
    epsilon = CONTRACT["current_IMEX"]["roundoff_epsilon_multiplier"]*np.finfo(float).eps
    limit = epsilon*steps*max(1, float(np.max(initial_h)))
    difference_h = actual-expected
    difference_eta = (actual+bed)-(expected+bed)
    result = {"cells_compared": actual.size, "limit_m": limit,
              "h_Linf_m": float(np.max(np.abs(difference_h))),
              "eta_Linf_m": float(np.max(np.abs(difference_eta))),
              "h_L1_m": float(np.mean(np.abs(difference_h))),
              "eta_L1_m": float(np.mean(np.abs(difference_eta)))}
    if max(result["h_Linf_m"], result["eta_Linf_m"]) > limit:
        raise AssertionError(f"current IMEX full-field mismatch: {result}")
    return result


def historical_comparison(actual, expected, bed):
    """Use explicitly approved non-roundoff limits only for the SSPRK2 trajectory."""
    result = {"normalized_L1_h": float(np.sum(np.abs(actual-expected))/np.sum(np.abs(expected))),
              "normalized_L1_eta": float(np.sum(np.abs((actual+bed)-(expected+bed)))/
                                         np.sum(np.abs(expected+bed))),
              "Linf_m": float(np.max(np.abs(actual-expected)))}
    limits = CONTRACT["historical_SSPRK2"]
    for field in result:
        if not np.isfinite(result[field]) or result[field] > limits["maximum_"+field]:
            raise AssertionError(f"historical SSPRK2 field tolerance failed: {result}")
    return result


def accepted_times(values, end):
    """Require finite positive contiguous accepted steps ending at the target time."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2 or len(values) == 0:
        raise AssertionError("invalid accepted-time stream layout")
    if not np.all(np.isfinite(values)) or np.any(values[:, 1] <= 0):
        raise AssertionError("nonfinite or nonpositive accepted dt")
    limit = 64*np.finfo(float).eps*max(1, abs(end))
    if abs(values[0, 0]) > limit or np.any(abs(values[1:, 0]-values[:-1].sum(axis=1)) > limit):
        raise AssertionError("accepted times are not a contiguous initial-value trajectory")
    if abs(values[-1].sum()-end) > limit:
        raise AssertionError("accepted trajectory did not reach comparison time")
    return values[:, 1]


def run(case_dir, mode, autonomous=False):
    """Generate both references from identical Q1 geometry; compare every final cell."""
    geometry = np.load(case_dir/"gate_h_geometry.npz")
    data = np.loadtxt(case_dir/"gateH_0001.q_2d")
    dx = float(geometry["dx"])
    bed = geometry["bed"]
    x, y = np.meshgrid(geometry["x"], geometry["y"])
    if data.shape != (bed.size, 6):
        raise AssertionError("invalid full-solver field layout")
    np.testing.assert_allclose(data[:, 0].reshape(bed.shape), x, rtol=0, atol=1e-10)
    np.testing.assert_allclose(data[:, 1].reshape(bed.shape), y, rtol=0, atol=1e-10)
    actual = data[:, 2].reshape(bed.shape)/1000
    flags = mode == "slope_curvature"
    settings = CONTRACT["current_IMEX"]
    initial = dynamic.equilibrium.liquid_state(geometry["h0"])
    time_evidence = json.loads((case_dir/"accepted_time_evidence.json").read_text())
    if not time_evidence["canonical_output_byte_identical"]:
        raise AssertionError("observed times are not certified against the actual output")
    for name, key in (("accepted_times.bin", "accepted_times_sha256"),
                      ("gateH_0001.q_2d", "canonical_output_sha256")):
        if hashlib.sha256((case_dir/name).read_bytes()).hexdigest() != time_evidence[key]:
            raise AssertionError("accepted-time evidence does not match the actual artifacts")
    dt = accepted_times(np.fromfile(case_dir/"accepted_times.bin").reshape(-1, 2), CONTRACT["time"])
    result = dynamic.reference(initial, geometry["bed_vertex"], dx, dx, dt,
        settings["n_RK"], flags, flags, rho_a=101300/(287.051*300),
        eps_sing=settings["eps_sing"])
    same = result["q"][..., 0]/1000
    failures = []
    try:
        current = compare(actual, same, bed, geometry["h0"], len(result["dt_sequence"]))
    except AssertionError as error:
        failures.append(str(error))
        current = {"status": "failed", "failure": str(error),
                   "h_Linf_m": float(np.max(np.abs(actual-same))),
                   "eta_Linf_m": float(np.max(np.abs((actual+bed)-(same+bed)))),
                   "cells_compared": actual.size}
    independent_time = {"status": "not_satisfied", "executed_in_this_run": False,
                        "closure_action": "D-N7-CFL", "N7_D_assertion": False,
                        "reason": "Autonomous adaptive trajectories are an explicitly separated unresolved diagnostic"}
    if autonomous:
        alternate = dynamic.reference(initial, geometry["bed_vertex"], dx, dx, None,
            settings["n_RK"], flags, flags, rho_a=101300/(287.051*300),
            eps_sing=settings["eps_sing"], time_policy={"dt0": settings["dt0"],
            "max_dt": settings["max_dt"], "end": CONTRACT["time"], "maximum_steps": 1000})
        try:
            assessment = compare(actual, alternate["q"][..., 0]/1000, bed,
                                 geometry["h0"], len(alternate["dt_sequence"]))
            independent_time.update(status="passed", assessment=assessment)
        except AssertionError as error:
            independent_time.update(status="failed", failure=str(error))
        independent_time.update(executed_in_this_run=True,
            dt_sequence=alternate["dt_sequence"], directional_bounds=alternate["bounds"])
    # Independently regenerated SSPRK2 fields, not the old port-audit scalar norms.
    settings = CONTRACT["historical_SSPRK2"]
    params = dynamic.core.SolverParams(rho_a=settings["rho_a"], eps_sing=settings["eps_sing"],
                                      slope_correction=flags, curvature_term=flags)
    bed_core = dynamic.core.continuous_q1_bed_geometry(geometry["bed_vertex"])
    q = dynamic.compact(initial); time = 0.0; prototype_steps = []
    while time < CONTRACT["time"]:
        q, diagnostic = dynamic.core.step_ssprk2_2d(q, bed_core, dx, dx, time,
            min(settings["dt_cap"], CONTRACT["time"]-time), params)
        prototype_steps.append(diagnostic)
        time += diagnostic["dt"]
    prototype = q[..., 0]/1000
    try:
        historical = historical_comparison(actual, prototype, bed)
    except AssertionError as error:
        failures.append(str(error))
        historical = {"status": "failed", "failure": str(error)}
    evidence = {"mode": mode, "dx": dx, "comparison_time": CONTRACT["time"],
        "contract_sha256": hashlib.sha256(CONTRACT_PATH.read_bytes()).hexdigest(),
        "reference_core_sha256": dynamic.spatial.CORE_HASHES,
        "current_IMEX": current, "current_IMEX_dt_sequence": result["dt_sequence"],
        "accepted_time_observer": time_evidence,
        "autonomous_CFL_diagnostic": independent_time,
        "current_IMEX_minimum_raw_h": result["minimum_raw_h"],
        "current_IMEX_directional_bounds": result["bounds"],
        "historical_SSPRK2": historical, "historical_SSPRK2_steps": prototype_steps,
        "reference_policy_equal": False, "old_port_norms_used_as_oracle": False,
        "status": "failed" if failures else "passed", "failures": failures}
    np.savez(case_dir/"reference_fields.npz", h_IMEX=same, eta_IMEX=same+bed,
             h_SSPRK2=prototype, eta_SSPRK2=prototype+bed)
    (case_dir/"field_evidence.json").write_text(json.dumps(evidence, indent=2, allow_nan=False)+"\n")
    if failures:
        raise AssertionError("; ".join(failures))
    print(f"PASS: Gate-H full fields {mode} dx={dx:g}, IMEX Linf={current['h_Linf_m']:.3e} m, "
          f"SSPRK2 Linf={historical['Linf_m']:.3e} m", flush=True)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("case_dir", type=Path)
    cli.add_argument("mode", choices=CONTRACT["modes"])
    cli.add_argument("--autonomous", action="store_true", help="Also record the unresolved autonomous CFL diagnostic")
    args = cli.parse_args()
    run(args.case_dir, args.mode, args.autonomous)
