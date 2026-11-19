#!/usr/bin/env python3
"""Localize Gate-H discrepancies using a certified accepted-step prefix.

Use the already built observational main (same numerical objects as the
production executable). Recreate only the inputs and stop at the exact requested
accepted-step prefix. Require matching observed steps and compare its conservative state
against independent KT-PYTHON forward-Euler updates.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np

from compare_fields import CONTRACT, accepted_times, common_input, kt


def replace_scalar(inp, name, value):
    pattern = rf"(?im)^(\s*{re.escape(name)}\s*=\s*)[^,\s]+"
    updated, count = re.subn(pattern, lambda match: match.group(1) + value, inp)
    if count != 1:
        raise AssertionError(f"expected one {name} setting, found {count}")
    return updated


def diagnose(observer, case_dir, mode, prefix_steps=1):
    case_dir = case_dir.resolve()
    obs = observer.resolve()
    geometry = np.load(case_dir / "gate_h_geometry.npz")
    dx = float(geometry["dx"])
    shape = geometry["bed"].shape
    times = np.fromfile(case_dir / "accepted_times.bin", dtype=np.float64).reshape(-1, 2)
    full_dt = accepted_times(times, CONTRACT["time"])
    if prefix_steps < 1 or prefix_steps > len(full_dt):
        raise ValueError("invalid accepted-step prefix length")
    target_dt = full_dt[:prefix_steps]
    stop_time = float(np.sum(target_dt))
    directory = case_dir / f"prefix-{prefix_steps}-step-diagnostic"
    directory.mkdir(exist_ok=True)
    for name in ("IMEX_SfloW2D.inp", "gateH_0000.q_2d"):
        shutil.copyfile(case_dir / name, directory / name)
    (directory / "DEM").mkdir(exist_ok=True)
    shutil.copyfile(case_dir / "DEM/topography_dem.asc",
                    directory / "DEM/topography_dem.asc")

    # Stop at the exact original accepted prefix endpoint; verify every replayed dt.
    time_value = f"{stop_time:.17e}".replace("e", "D")
    inp = (directory / "IMEX_SfloW2D.inp").read_text()
    inp = replace_scalar(inp, "T_END", time_value)
    inp = replace_scalar(inp, "DT_OUTPUT", time_value)
    (directory / "IMEX_SfloW2D.inp").write_text(inp)
    with (directory / "run.log").open("w") as log:
        subprocess.run([str(obs)], cwd=directory, env=os.environ.copy(),
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=600)

    step_file = directory / "accepted_times.bin"
    if not step_file.exists():
        raise AssertionError("first-step observer did not record any accepted times")
    steps = np.fromfile(step_file, dtype=np.float64).reshape(-1, 2)
    if steps.shape != (prefix_steps, 2):
        raise AssertionError(f"expected {prefix_steps} steps, observed {steps.tolist()}")
    if not np.allclose(steps[:, 1], target_dt, rtol=0, atol=16*np.finfo(float).eps):
        raise AssertionError(f"replayed accepted dt differs: {steps[:, 1]} != {target_dt}")
    expected_t = np.r_[0.0, np.cumsum(target_dt)[:-1]]
    if not np.allclose(steps[:, 0], expected_t, rtol=0, atol=16*np.finfo(float).eps):
        raise AssertionError(f"replayed accepted start times differ: {steps[:, 0]} != {expected_t}")

    snapshots = sorted((directory / "OUTPUT/gateH/snapshots").glob("*/gateH_0001.q_2d"))
    if len(snapshots) != 1:
        raise AssertionError(f"expected exactly one first-step snapshot: {snapshots}")
    data = np.loadtxt(snapshots[0])
    if data.shape != (int(np.prod(shape)), 6):
        raise AssertionError(f"invalid replayed snapshot shape: {data.shape}")
    actual = data[:, 2:5].reshape(shape + (3,))

    initial, vertices = common_input(case_dir, shape, dx)
    flags = mode == "slope_curvature"
    settings = CONTRACT["current_IMEX"]
    params = kt.Parameters(rho_a=101300/(287.051*300),
                           eps_sing=settings["eps_sing"],
                           reconstruction=settings["KT_TRACE_RECONSTRUCTION"],
                           slope_correction=flags, curvature_term=flags)
    residual = kt.rhs(initial, vertices, dx, dx, params)
    reference, _, _ = kt.advance(initial, vertices, dx, dx, target_dt, params,
                                  return_diagnostics=True)
    expected = reference[..., :3]
    differences = actual - expected
    scales = [1000.0, 1000.0*np.sqrt(2*params.g), 1000.0*np.sqrt(2*params.g)]
    normalized = np.max(np.abs(differences), axis=(0, 1)) / scales

    # The full-prefix mean rate is diagnostic; it equals the spatial residual
    # only for a single step. For two or more steps, compare conserved fields.
    observed_residual = (actual-initial[..., :3])/stop_time
    residual_difference = observed_residual - residual[..., :3]
    worst = np.unravel_index(np.argmax(np.abs(differences[..., 0])), shape)
    j, i = (int(index) for index in worst)
    patch = (slice(max(0,j-1), min(shape[0],j+2)),
             slice(max(0,i-1), min(shape[1],i+2)))
    report = {
        "mode": mode, "dx": dx, "prefix_steps": prefix_steps,
        "accepted_prefix_dt": target_dt.tolist(), "stop_time": stop_time,
        "original_full_run_steps": int(len(full_dt)),
        "observed_prefix_dt": steps[:, 1].tolist(),
        "max_absolute_conservative_error": np.max(np.abs(differences), axis=(0,1)).tolist(),
        "max_scaled_conservative_error": normalized.tolist(),
        "max_absolute_inferred_residual_error": np.max(
            np.abs(residual_difference), axis=(0,1)).tolist(),
        "worst_mass_cell_ij": [i,j],
        "worst_mass_cell_initial_q": initial[j,i,:3].tolist(),
        "worst_mass_cell_fortran_q": actual[j,i,:].tolist(),
        "worst_mass_cell_python_q": expected[j,i,:].tolist(),
        "worst_mass_cell_inferred_fortran_residual": observed_residual[j,i,:].tolist(),
        "worst_mass_cell_python_residual": residual[j,i,:3].tolist(),
        "mass_error_3x3_kg_m2": differences[patch][...,0].tolist(),
        "momentum_x_error_3x3": differences[patch][...,1].tolist(),
        "momentum_y_error_3x3": differences[patch][...,2].tolist(),
        "meaning": ("A matching prefix followed by a mismatching prefix isolates "
                    "the first divergent time interval. Mean rates over several "
                    "steps must not be compared to the initial spatial residual."),
    }
    out = case_dir / f"prefix_{prefix_steps}_step_evidence.json"
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"PREFIX-{prefix_steps}-STEP {mode} dx={dx:g} stop={stop_time:.17g} "
          f"max-scaled-error={normalized.tolist()} worst cell (i,j)=({i},{j})",flush=True)
    print(f"PREFIX-{prefix_steps}-STEP mean residual Fortran={observed_residual[j,i,:].tolist()} "
          f"Python={residual[j,i,:3].tolist()}",flush=True)
    print(f"PREFIX-{prefix_steps}-STEP evidence={out}",flush=True)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("observer", type=Path)
    cli.add_argument("case_dir", type=Path)
    cli.add_argument("mode", choices=CONTRACT["modes"])
    cli.add_argument("--steps", type=int, default=1,
                     help="length of replayed accepted-step prefix")
    args = cli.parse_args()
    diagnose(args.observer, args.case_dir, args.mode, args.steps)
