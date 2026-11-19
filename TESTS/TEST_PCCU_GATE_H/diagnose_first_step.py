#!/usr/bin/env python3
"""Localize Gate-H discrepancies at the first accepted step.

Use the already built observational main (same numerical objects as the
production executable). Recreate only the inputs and stop at the exact first
accepted dt. Require one observed step and compare its conservative state
against an independent KT-PYTHON forward-Euler update.
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


def diagnose(observer, case_dir, mode):
    case_dir = case_dir.resolve()
    obs = observer.resolve()
    geometry = np.load(case_dir / "gate_h_geometry.npz")
    dx = float(geometry["dx"])
    shape = geometry["bed"].shape
    times = np.fromfile(case_dir / "accepted_times.bin", dtype=np.float64).reshape(-1, 2)
    full_dt = accepted_times(times, CONTRACT["time"])
    first_dt = float(full_dt[0])
    directory = case_dir / "first-step-diagnostic"
    directory.mkdir(exist_ok=True)
    for name in ("IMEX_SfloW2D.inp", "gateH_0000.q_2d"):
        shutil.copyfile(case_dir / name, directory / name)
    (directory / "DEM").mkdir(exist_ok=True)
    shutil.copyfile(case_dir / "DEM/topography_dem.asc",
                    directory / "DEM/topography_dem.asc")

    # The final time and output event coincide with the original first step.
    # The observer records the actual step taken so it can be checked exactly.
    time_value = f"{first_dt:.17e}".replace("e", "D")
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
    if steps.shape != (1, 2) or not np.isclose(steps[0, 0], 0.0, rtol=0, atol=1e-15):
        raise AssertionError(f"expected exactly one initial step, observed {steps.tolist()}")
    if not np.isclose(steps[0, 1], first_dt, rtol=0, atol=16*np.finfo(float).eps):
        raise AssertionError(f"first-step replay dt={steps[0, 1]} differs from full-run dt={first_dt}")

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
    reference, _, _ = kt.advance(initial, vertices, dx, dx, [first_dt], params,
                                  return_diagnostics=True)
    expected = reference[..., :3]
    differences = actual - expected
    scales = [1000.0, 1000.0*np.sqrt(2*params.g), 1000.0*np.sqrt(2*params.g)]
    normalized = np.max(np.abs(differences), axis=(0, 1)) / scales

    # Residual inferred from the observed one-step conservative increment.
    # For n_RK=2 on this source-free case, the spatial part is forward Euler.
    observed_residual = (actual-initial[..., :3])/first_dt
    residual_difference = observed_residual - residual[..., :3]
    worst = np.unravel_index(np.argmax(np.abs(differences[..., 0])), shape)
    j, i = (int(index) for index in worst)
    patch = (slice(max(0,j-1), min(shape[0],j+2)),
             slice(max(0,i-1), min(shape[1],i+2)))
    report = {
        "mode": mode, "dx": dx, "first_dt": first_dt,
        "original_full_run_steps": int(len(full_dt)),
        "first_step_observed_dt": float(steps[0, 1]),
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
        "meaning": ("Initial-step difference isolates the first explicit spatial/source "
                    "update; a near-zero first-step difference would instead implicate "
                    "subsequent state updates or time evolution."),
    }
    out = case_dir / "first_step_evidence.json"
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"FIRST-STEP {mode} dx={dx:g} dt={first_dt:.17g} "
          f"max-scaled-error={normalized.tolist()} worst cell (i,j)=({i},{j})",flush=True)
    print(f"FIRST-STEP local residual Fortran={observed_residual[j,i,:].tolist()} "
          f"Python={residual[j,i,:3].tolist()}",flush=True)
    print(f"FIRST-STEP evidence={out}",flush=True)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("observer", type=Path)
    cli.add_argument("case_dir", type=Path)
    cli.add_argument("mode", choices=CONTRACT["modes"])
    args = cli.parse_args()
    diagnose(args.observer, args.case_dir, args.mode)
