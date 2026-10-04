#!/usr/bin/env python3
"""Check Gate-H invariants and compact frozen-reference field fingerprints."""

from pathlib import Path
import argparse
import csv
import numpy as np


def load_reference(path, mode, dx):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    matches = [
        row
        for row in rows
        if row["mode"] == mode and np.isclose(float(row["dx"]), dx)
    ]
    if len(matches) != 1:
        raise RuntimeError(f"missing unique reference row for mode={mode}, dx={dx}")
    return {
        key: (value if key == "mode" else float(value))
        for key, value in matches[0].items()
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("mode", choices=("g1", "slope_curvature"))
    parser.add_argument("reference", type=Path)
    args = parser.parse_args()

    geometry = np.load(args.case_dir / "gate_h_geometry.npz")
    state = np.loadtxt(args.case_dir / "gateH_0001.q_2d")
    dx = float(geometry["dx"])
    x_grid = geometry["x"]
    y_grid = geometry["y"]
    nx = x_grid.size
    ny = y_grid.size
    if state.shape[0] != nx * ny:
        raise SystemExit(f"unexpected output size: {state.shape[0]} != {nx * ny}")

    x = state[:, 0]
    y = state[:, 1]
    mass = state[:, 2]
    h = (mass / 1000.0).reshape(ny, nx)
    X, Y = np.meshgrid(x_grid, y_grid)
    initial_mass = 1000.0 * float(np.sum(geometry["h0"])) * dx * dx
    final_mass = float(np.sum(mass)) * dx * dx
    mass_relative_error = (final_mass - initial_mass) / initial_mass

    x_left = float(geometry["x_outer_left"])
    x_right = float(geometry["x_outer_right"])
    y_bottom = float(geometry["y_outer_bottom"])
    y_top = float(geometry["y_outer_top"])
    region_masks = {
        "uphill": x < x_left,
        "downstream": x >= x_right,
        "bottom": (x >= x_left) & (x < x_right) & (y < y_bottom),
        "top": (x >= x_left) & (x < x_right) & (y >= y_top),
    }
    fractions = {
        name: float(np.sum(mass[mask]) * dx * dx / initial_mass)
        for name, mask in region_masks.items()
    }
    fractions["lateral"] = fractions["bottom"] + fractions["top"]

    h_sum = max(float(np.sum(h)), 1.0e-300)
    metrics = {
        "hmax": float(np.max(h)),
        "h_l2": float(np.sqrt(np.mean(h * h))),
        "h_xmean": float(np.sum(h * X) / h_sum),
        "h_ymean": float(np.sum(h * Y) / h_sum),
        "h_x2mean": float(np.sum(h * X * X) / h_sum),
        "h_y2mean": float(np.sum(h * Y * Y) / h_sum),
        "tvx": float(np.sum(np.abs(np.diff(h, axis=1))) / h_sum),
        "tvy": float(np.sum(np.abs(np.diff(h, axis=0))) / h_sum),
    }
    reference = load_reference(args.reference, args.mode, dx)

    failures = []
    if float(np.min(mass)) < -1.0e-10:
        failures.append(f"negative mass {np.min(mass):.3e}")
    if abs(mass_relative_error) > 5.0e-10:
        failures.append(f"mass relative error {mass_relative_error:.3e}")
    if abs(fractions["uphill"]) > 1.0e-14:
        failures.append(f"uphill fraction {fractions['uphill']:.3e}")

    partition_tolerance = {
        "downstream": 1.0e-5 if args.mode == "g1" else 4.0e-6,
        "lateral": 2.0e-6 if args.mode == "g1" else 1.0e-7,
    }
    for name, tolerance in partition_tolerance.items():
        error = abs(fractions[name] - reference[name])
        if error > tolerance:
            failures.append(
                f"{name} reference error {error:.3e} > {tolerance:.3e}"
            )

    # These aggregate fingerprints make the regression sensitive to the full
    # thickness field without storing three large binary reference arrays.
    tolerances = {
        "hmax": 3.0e-3,
        "h_l2": 1.0e-3,
        "h_xmean": 5.0e-3,
        "h_ymean": 2.0e-5,
        "h_x2mean": 8.0e-2,
        "h_y2mean": 2.0e-2,
        "tvx": 3.0e-3,
        "tvy": 3.0e-3,
    }
    for name, tolerance in tolerances.items():
        error = abs(metrics[name] - reference[name])
        if error > tolerance:
            failures.append(
                f"{name} reference error {error:.3e} > {tolerance:.3e}"
            )

    # With G=1 the compact reference retains y symmetry to much tighter than
    # the field-comparison tolerances.  The slope LS stencil is intentionally
    # not checked here because the coarse compact domain is phase-asymmetric
    # with respect to its copied boundary strip.
    if args.mode == "g1":
        symmetry_error = abs(fractions["bottom"] - fractions["top"])
        if symmetry_error > 5.0e-10:
            failures.append(f"lateral symmetry error {symmetry_error:.3e}")

    print(
        f"mode={args.mode} dx={dx:g} mass_error={mass_relative_error:.3e} "
        f"uphill={fractions['uphill']:.3e} lateral={fractions['lateral']:.3e} "
        f"downstream={fractions['downstream']:.3e} hmax={metrics['hmax']:.9e}"
    )
    if failures:
        raise SystemExit("; ".join(failures))


if __name__ == "__main__":
    main()
