"""Compare all Gate-H cells against the current independent KT-PYTHON scheme."""
import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REFERENCE = HERE.parent / "TEST_SPATIAL_OPERATOR/reference"
sys.path.insert(0, str(REFERENCE))
import kt_python_2d_core as kt
CONTRACT_PATH = HERE / "field_contract_v4.json"
CONTRACT = json.loads(CONTRACT_PATH.read_text())
for filename, digest in CONTRACT["reference_core_sha256"].items():
    if hashlib.sha256((REFERENCE / filename).read_bytes()).hexdigest() != digest:
        raise AssertionError(f"reference implementation checksum mismatch: {filename}")


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


def common_input(case_dir, shape, dx):
    """Independently import the exact ASCII initial state and sampled DEM.

    Reading the emitted input bytes avoids comparing unrounded generator
    arrays with the rounded data actually read by the production executable.
    """
    initial_data = np.loadtxt(case_dir/"gateH_0000.q_2d")
    state = np.zeros(shape+(4,))
    state[..., :3] = initial_data[:,2:5].reshape(shape+(3,))
    with (case_dir/"topography_dem.asc").open() as stream:
        header = dict((key.lower(),float(value)) for key,value in
                      (stream.readline().split() for _ in range(6)))
        raster = np.flipud(np.loadtxt(stream))
    inp = (case_dir/"IMEX_SfloW2D.inp").read_text()
    def value(name):
        return float(re.search(rf"(?im)^\s*{name}\s*=\s*([^,\s]+)",inp)[1].replace('D','E').replace('d','e'))
    x = value('X0')+np.arange(shape[1]+1)*dx
    y = value('Y0')+np.arange(shape[0]+1)*dx
    xd = header['xllcorner']+(np.arange(raster.shape[1])+0.5)*header['cellsize']
    yd = header['yllcorner']+(np.arange(raster.shape[0])+0.5)*header['cellsize']
    ix = np.clip(np.floor((x-xd[0])/(xd[1]-xd[0])).astype(int),0,len(xd)-2)
    iy = np.clip(np.floor((y-yd[0])/(yd[1]-yd[0])).astype(int),0,len(yd)-2)
    ax = (xd[ix+1]-x)/(xd[ix+1]-xd[ix])
    ay = (yd[iy+1]-y)/(yd[iy+1]-yd[iy])
    vertices = ax[None,:]*(ay[:,None]*raster[iy[:,None],ix[None,:]]+
        (1-ay[:,None])*raster[iy[:,None]+1,ix[None,:]])+(1-ax[None,:])*(
        ay[:,None]*raster[iy[:,None],ix[None,:]+1]+
        (1-ay[:,None])*raster[iy[:,None]+1,ix[None,:]+1])
    return state,vertices


def run(case_dir, mode, autonomous=False):
    """Generate the current KT reference from common inputs; compare every cell."""
    if autonomous:
        raise ValueError("autonomous CFL comparison is not implemented for the new KT core")
    geometry = np.load(case_dir/"gate_h_geometry.npz")
    snapshot = next((case_dir/"OUTPUT/gateH/snapshots").glob("*/gateH_0001.q_2d"))
    data = np.loadtxt(snapshot)
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
    initial, vertices = common_input(case_dir,bed.shape,dx)
    time_evidence = json.loads((case_dir/"accepted_time_evidence.json").read_text())
    if not time_evidence["canonical_output_byte_identical"]:
        raise AssertionError("observed times are not certified against the actual output")
    for path, key in ((case_dir/"accepted_times.bin", "accepted_times_sha256"),
                      (snapshot, "canonical_output_sha256")):
        if hashlib.sha256(path.read_bytes()).hexdigest() != time_evidence[key]:
            raise AssertionError("accepted-time evidence does not match the actual artifacts")
    dt = accepted_times(np.fromfile(case_dir/"accepted_times.bin").reshape(-1, 2), CONTRACT["time"])
    if settings["n_RK"] != 2:
        raise AssertionError("this comparison requires n_RK=2")
    params = kt.Parameters(rho_a=101300/(287.051*300), eps_sing=settings["eps_sing"],
                           reconstruction=settings["KT_TRACE_RECONSTRUCTION"],
                           slope_correction=flags, curvature_term=flags)
    reference, minimum, budget = kt.advance(initial, vertices, dx, dx, dt, params,
                                               return_diagnostics=True)
    same = reference[..., 0]/1000
    failures = []
    try:
        current = compare(actual, same, bed, geometry["h0"], len(dt))
        # Check momenta too: equal thickness alone must not hide a path/flux error.
        actual_q = data[:, 2:5].reshape(bed.shape+(3,))
        # Scale each equation by the independent initial hydrostatic units,
        # also for almost stationary/near-dry momentum (N7 convention).
        mass_scale = max(1,float(initial[...,0].max()))
        momentum_scale = mass_scale*np.sqrt(params.g*max(1,float(initial[...,0].max())/params.rho_c))
        scales = np.array([mass_scale,momentum_scale,momentum_scale])
        limit = settings["roundoff_epsilon_multiplier"]*np.finfo(float).eps*len(dt)
        errors = np.max(abs(actual_q-reference[..., :3]), axis=(0, 1))/scales
        if np.any(errors > limit):
            raise AssertionError(f"KT conservative field mismatch: {errors}, limit={limit}")
        current["conservative_scaled_Linf"] = errors.tolist()
        current["conservative_scaled_limit"] = limit
    except AssertionError as error:
        failures.append(str(error))
        current = {"status": "failed", "failure": str(error),
                   "h_Linf_m": float(np.max(np.abs(actual-same))),
                   "cells_compared": actual.size}
    evidence = {"mode": mode, "dx": dx, "comparison_time": CONTRACT["time"],
        "contract_sha256": hashlib.sha256(CONTRACT_PATH.read_bytes()).hexdigest(),
        "reference_core_sha256": CONTRACT["reference_core_sha256"],
        "scheme": "KT-PYTHON", "trace_reconstruction": params.reconstruction,
        "n_RK": 2, "boundary_budget": budget, "current_IMEX": current, "current_IMEX_dt_sequence": dt.tolist(),
        "accepted_time_observer": time_evidence, "current_IMEX_minimum_raw_h": minimum,
        "historical_SSPRK2": {"acceptance_oracle": False,
            "reason": "different spatial scheme; retained historical contract is not applied"},
        "autonomous_CFL_diagnostic": {"executed_in_this_run": False,
            "reason": "this comparison uses the certified actual accepted time sequence"},
        "status": "failed" if failures else "passed", "failures": failures}
    np.savez(case_dir/"reference_fields.npz", h_IMEX=same, eta_IMEX=same+bed, q_IMEX=reference)
    evidence["reference_fields_sha256"] = hashlib.sha256((case_dir/"reference_fields.npz").read_bytes()).hexdigest()
    (case_dir/"field_evidence.json").write_text(json.dumps(evidence, indent=2, allow_nan=False)+"\n")
    if failures:
        raise AssertionError("; ".join(failures))
    print(f"PASS: KT-PYTHON full fields {mode} dx={dx:g}, h Linf={current['h_Linf_m']:.3e} m", flush=True)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("case_dir", type=Path)
    cli.add_argument("mode", choices=CONTRACT["modes"])
    cli.add_argument("--autonomous", action="store_true", help="Retained legacy option; reports that KT autonomous CFL comparison is unavailable")
    args = cli.parse_args()
    run(args.case_dir, args.mode, args.autonomous)
