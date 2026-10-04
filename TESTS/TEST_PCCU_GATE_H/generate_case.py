#!/usr/bin/env python3
"""Generate the compact four-sided one-cell Gate-H excavation."""

from pathlib import Path
import argparse
import numpy as np


def ramp_indicator(coordinate, left_center, right_center, width):
    """Continuous one-cell ramp with walls centred at the supplied points."""
    left = np.clip(
        (coordinate - (left_center - 0.5 * width)) / width, 0.0, 1.0
    )
    right = np.clip(
        ((right_center + 0.5 * width) - coordinate) / width, 0.0, 1.0
    )
    return np.minimum(left, right)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dx", type=float)
    parser.add_argument("mode", choices=("g1", "slope_curvature"))
    parser.add_argument("outdir", type=Path)
    parser.add_argument("template", type=Path)
    args = parser.parse_args()

    dx = args.dx
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    x_min, x_max = 10.8, 21.2
    y_min, y_max = -4.2, 4.2
    nx = round((x_max - x_min) / dx)
    ny = round((y_max - y_min) / dx)
    if not np.isclose(x_min + nx * dx, x_max):
        raise ValueError("dx does not divide the Gate-H x domain")
    if not np.isclose(y_min + ny * dx, y_max):
        raise ValueError("dx does not divide the Gate-H y domain")

    x_center = x_min + (np.arange(nx) + 0.5) * dx
    y_center = y_min + (np.arange(ny) + 0.5) * dx
    x_vertex = x_min + np.arange(nx + 1) * dx
    y_vertex = y_min + np.arange(ny + 1) * dx
    Xc, Yc = np.meshgrid(x_center, y_center)
    Xv, Yv = np.meshgrid(x_vertex, y_vertex)

    # The nominal 8 m x 6 m rectangle is translated by half a cell.  Its
    # one-cell ramps therefore occupy [x1,x1+dx], [x2,x2+dx] and analogously
    # in y, exactly as in the frozen Python Gate-H benchmark.
    x1, x2 = 12.0, 20.0
    y_halfwidth = 3.0
    depth = 2.0
    slope = np.tan(np.deg2rad(10.0))
    z_reference = 10.0

    Dx = ramp_indicator(Xv, x1 + 0.5 * dx, x2 + 0.5 * dx, dx)
    Dy = ramp_indicator(
        Yv, -y_halfwidth + 0.5 * dx, y_halfwidth + 0.5 * dx, dx
    )
    depth_vertex = depth * Dx * Dy
    bed_vertex = z_reference - slope * Xv - depth_vertex
    bed_center = 0.25 * (
        bed_vertex[:-1, :-1]
        + bed_vertex[1:, :-1]
        + bed_vertex[:-1, 1:]
        + bed_vertex[1:, 1:]
    )
    thickness = 0.25 * (
        depth_vertex[:-1, :-1]
        + depth_vertex[1:, :-1]
        + depth_vertex[:-1, 1:]
        + depth_vertex[1:, 1:]
    )

    # Add one DEM sample outside each computational boundary.  The inner DEM
    # samples then coincide with the authoritative computational vertices.
    x_dem = x_min - dx + np.arange(nx + 3) * dx
    y_dem = y_min - dx + np.arange(ny + 3) * dx
    Xd, Yd = np.meshgrid(x_dem, y_dem)
    Dx_dem = ramp_indicator(Xd, x1 + 0.5 * dx, x2 + 0.5 * dx, dx)
    Dy_dem = ramp_indicator(
        Yd, -y_halfwidth + 0.5 * dx, y_halfwidth + 0.5 * dx, dx
    )
    bed_dem = z_reference - slope * Xd - depth * Dx_dem * Dy_dem
    header = (
        f"ncols {nx + 3}\n"
        f"nrows {ny + 3}\n"
        f"xllcorner {x_min - 1.5 * dx}\n"
        f"yllcorner {y_min - 1.5 * dx}\n"
        f"cellsize {dx}\n"
        "NODATA_value -9999\n"
    )
    np.savetxt(
        outdir / "topography_dem.asc",
        np.flipud(bed_dem),
        header=header,
        comments="",
        fmt="%.15e",
    )

    rho_liquid = 1000.0
    specific_heat = 4200.0
    temperature = 300.0
    with (outdir / "gateH_0000.q_2d").open("w") as stream:
        for j in range(ny):
            state = np.column_stack(
                (
                    Xc[j],
                    Yc[j],
                    rho_liquid * thickness[j],
                    np.zeros(nx),
                    np.zeros(nx),
                    rho_liquid * thickness[j] * specific_heat * temperature,
                )
            )
            np.savetxt(stream, state, fmt="%19.12e")
            stream.write(" \n")

    inp = args.template.read_text()
    replacements = {
        "runname": "gateH",
        "restartfile": "gateH_0000.q_2d",
        "x_min": str(x_min),
        "y_min": str(y_min),
        "nx_cells": str(nx),
        "ny_cells": str(ny),
        "dx": repr(dx),
        "T_END=3.0D0": "T_END=1.0D-1",
        "DT_OUTPUT=5.0D-2": "DT_OUTPUT=1.0D-1",
        "OUTPUT_CONS_FLAG=F": "OUTPUT_CONS_FLAG=T",
        "OUTPUT_NETCDF_FLAG=T": "OUTPUT_NETCDF_FLAG=F",
    }
    if args.mode == "slope_curvature":
        replacements.update(
            {
                "SLOPE_CORRECTION_FLAG=F": "SLOPE_CORRECTION_FLAG=T",
                "CURVATURE_TERM_FLAG=F": "CURVATURE_TERM_FLAG=T",
            }
        )
    for old, new in replacements.items():
        if old not in inp:
            raise RuntimeError(f"missing template token: {old}")
        inp = inp.replace(old, new)
    (outdir / "IMEX_SfloW2D.inp").write_text(inp)

    np.savez(
        outdir / "gate_h_geometry.npz",
        x=x_center,
        y=y_center,
        bed=bed_center,
        h0=thickness,
        dx=dx,
        x_outer_left=x1,
        x_outer_right=x2 + dx,
        y_outer_bottom=-y_halfwidth,
        y_outer_top=y_halfwidth + dx,
    )

    mass = rho_liquid * np.sum(thickness) * dx * dx
    print(f"mode={args.mode} dx={dx:g} nx={nx} ny={ny} initial_mass={mass:.15e}")


if __name__ == "__main__":
    main()
