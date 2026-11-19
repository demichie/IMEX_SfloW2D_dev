#!/usr/bin/env python3
"""Inspect mass and wet-front speeds in the scalar-field NetCDF output.

This is an observational diagnostic, not a conservation acceptance test.
Open boundaries, injection, deposition and other sources can change inventory.
"""
import argparse
import csv
from pathlib import Path

import numpy as np


def open_fields(path):
    """Return a dataset and field reader; prefer netCDF4, then h5py."""
    try:
        from netCDF4 import Dataset

        dataset = Dataset(str(path), "r")
        def read(name, step=None):
            field = dataset.variables[name]
            return np.asarray(field[:] if step is None else field[step], dtype=float)
        return dataset, set(dataset.variables), read
    except ImportError:
        import h5py

        dataset = h5py.File(path, "r")
        def read(name, step=None):
            field = dataset[name]
            return np.asarray(field[...] if step is None else field[step], dtype=float)
        return dataset, set(dataset.keys()), read


def inspect(path, threshold):
    """Compute snapshot inventories and high-speed cell diagnostics."""
    dataset, names, read = open_fields(path)
    try:
        required = {"x", "y", "time", "h", "u", "v", "mass"}
        missing = required - names
        if missing:
            raise ValueError("Missing physical NetCDF fields: " + ", ".join(sorted(missing)))
        x, y, times = read("x"), read("y"), read("time")
        if len(x) < 2 or len(y) < 2:
            raise ValueError("At least two cells per horizontal axis are required")
        dx, dy = np.diff(x), np.diff(y)
        if not np.allclose(dx, dx[0]) or not np.allclose(dy, dy[0]):
            raise ValueError("This diagnostic expects a uniform Cartesian grid")
        area = float(abs(dx[0] * dy[0]))
        rows = []
        for step, time in enumerate(times):
            h, u, v, mass = (read(name, step) for name in ("h", "u", "v", "mass"))
            if h.shape != (len(y), len(x)):
                raise ValueError("Unexpected field shape: " + str(h.shape))
            speed = np.hypot(u, v)
            significant = h > threshold
            if np.any(significant):
                valid_speed = np.where(significant, speed, -1.0)
                index = np.unravel_index(int(np.argmax(valid_speed)), h.shape)
                max_speed = float(speed[index])
                depth_at_max = float(h[index])
                x_at_max = float(x[index[1]])
                y_at_max = float(y[index[0]])
            else:
                max_speed = depth_at_max = x_at_max = y_at_max = float("nan")
            rows.append({
                "step": step,
                "time_s": float(time),
                "mass_kg": float(area * np.sum(mass)),
                "volume_m3": float(area * np.sum(h)),
                "max_speed_m_s": max_speed,
                "depth_at_max_speed_m": depth_at_max,
                "x_at_max_speed_m": x_at_max,
                "y_at_max_speed_m": y_at_max,
                "min_positive_depth_m": float(np.min(h[h > 0])) if np.any(h > 0) else float("nan"),
                "wet_cells": int(np.count_nonzero(significant)),
                "max_pore_pressure_Pa": float(np.max(read("pore_pressure", step)))
                if "pore_pressure" in names else float("nan"),
            })
        return rows
    finally:
        dataset.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("netcdf", type=Path, help="Scalar-field NetCDF output")
    parser.add_argument("--wet-depth", type=float, default=1e-6,
                        help="Minimum depth for velocity diagnostics (m)")
    parser.add_argument("--csv", type=Path, default=None,
                        help="Output CSV (default: alongside input with .diagnostics.csv suffix)")
    args = parser.parse_args()
    if args.wet_depth < 0:
        parser.error("--wet-depth must be nonnegative")
    rows = inspect(args.netcdf, args.wet_depth)
    if not rows:
        parser.error("NetCDF contains no time records")
    destination = args.csv or args.netcdf.with_suffix(".diagnostics.csv")
    with destination.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    mass0 = rows[0]["mass_kg"]
    mass1 = rows[-1]["mass_kg"]
    relative_change = (mass1 - mass0) / mass0 if mass0 != 0 else float("nan")
    fastest = max(rows, key=lambda item: item["max_speed_m_s"]
                  if np.isfinite(item["max_speed_m_s"]) else -1.0)
    print("Records:", len(rows), "last_time_s:", rows[-1]["time_s"])
    print("Mass at first/last output (kg):", mass0, mass1)
    print("Relative mass change (NOT a conservation residual):", relative_change)
    print("Maximum resolved speed (m/s):", fastest["max_speed_m_s"],
          "at t=", fastest["time_s"], "h=", fastest["depth_at_max_speed_m"])
    print("CSV:", destination)
    print("Conservation requires accounting for boundary fluxes and physical mass sources.")


if __name__ == "__main__":
    main()
