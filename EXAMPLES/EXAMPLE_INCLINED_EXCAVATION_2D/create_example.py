#!/usr/bin/env python3
"""Create a 2-D inclined-excavation wet/dry test for IMEX_SfloW2D_dev.

The unexcavated bed is a plane descending in +x:
    B0(x) = z_ref - tan(slope_angle) * x.

Inside a rectangular patch, a constant thickness d is excavated:
    B = B0 - d.

The excavation is filled with pure water to h=d, so the initial free surface is
    eta = B + h = B0
inside the excavation. Thus eta is the same constant-slope plane as the
original unexcavated topography.

Expected physics:
  * initial acceleration is in +x (downslope);
  * no physical acceleration in y;
  * initially no water should cross the upslope/west, north, or south vertical
    excavation walls;
  * water should move internally toward the east/downhill wall and eventually
    overtop it.

The DEM contains a one-cell ghost ring around the computational domain.  Its
interior raster sample centers coincide exactly with the computational cell
centers, so the excavation step remains aligned with cell interfaces while the
DEM still extends far enough to satisfy the domain-coverage checks in
IMEX_SfloW2D.
"""

from pathlib import Path
import argparse
import numpy as np


def write_esri_cell_centered(path, z_comp, x_min, y_min, dx, slope, z_ref,
                             x1, x2, y1, y2, depth):
    """Write a padded ESRI ASCII DEM aligned with computational cell centers.

    IMEX_SfloW2D requires the DEM to cover the computational *boundaries*, not
    only the computational cell centers.  A DEM whose first raster center is at
    x_min + dx/2 therefore fails the check

        x0 >= xllcorner + 0.5*cellsize.

    We solve this by adding one DEM sample center on every side.  The DEM sample
    centers are then

        x_min-dx/2, x_min+dx/2, ..., x_max-dx/2, x_max+dx/2,

    and similarly in y.  The interior DEM centers coincide exactly with the
    computational centers, so no interpolation smears the excavation step.
    """
    ny, nx = z_comp.shape

    x_dem = x_min - 0.5 * dx + np.arange(nx + 2) * dx
    y_dem = y_min - 0.5 * dx + np.arange(ny + 2) * dx
    Xd, Yd = np.meshgrid(x_dem, y_dem)

    B0_dem = z_ref - slope * Xd
    mask_dem = (Xd >= x1) & (Xd < x2) & (Yd >= y1) & (Yd < y2)
    z_dem = B0_dem.copy()
    z_dem[mask_dem] -= depth

    # Sanity check: the interior DEM samples must reproduce the computational
    # topography exactly.
    if not np.allclose(z_dem[1:-1, 1:-1], z_comp, rtol=0.0, atol=1.0e-12):
        raise RuntimeError("Padded DEM interior is not aligned with computational cells")

    header = f"ncols     {nx + 2}\n"
    header += f"nrows    {ny + 2}\n"
    # ESRI xllcorner/yllcorner refer to the lower-left raster EDGE.  Since the
    # first DEM center is x_min-dx/2, the lower-left edge is x_min-dx.
    header += f"xllcorner {x_min - dx}\n"
    header += f"yllcorner {y_min - dx}\n"
    header += f"cellsize {dx}\n"
    header += "NODATA_value -9999\n"

    # ESRI ASCII rows are north -> south; the solver reverses them on input.
    np.savetxt(path, np.flipud(z_dem), header=header, fmt="%1.12f", comments="")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("nx", nargs="?", type=int, default=120,
                        help="number of x cells (default: 120)")
    parser.add_argument("--plot", action="store_true",
                        help="save diagnostic PNGs")
    parser.add_argument("--slope-deg", type=float, default=10.0,
                        help="bed slope angle in degrees, downhill toward +x")
    parser.add_argument("--depth", type=float, default=2.0,
                        help="constant excavation/fill depth [m]")
    args = parser.parse_args()

    # Fluid properties: pure water, no particles.
    rho_l = 1000.0
    sp_heat_l = 4200.0
    temperature = 300.0

    # Domain follows the scale of EXAMPLE_2D.
    x_min, x_max = 0.0, 30.0
    y_target_min, y_target_max = -10.0, 10.0

    nx_cells = args.nx
    dx = (x_max - x_min) / nx_cells

    # Square computational cells, as in EXAMPLE_2D.
    ny_cells = int(round((y_target_max - y_target_min) / dx))
    if ny_cells < 4:
        raise ValueError("Need at least 4 cells in y")
    y_min = -0.5 * ny_cells * dx
    y_max = +0.5 * ny_cells * dx

    x_cent = x_min + (np.arange(nx_cells) + 0.5) * dx
    y_cent = y_min + (np.arange(ny_cells) + 0.5) * dx
    Xc, Yc = np.meshgrid(x_cent, y_cent)

    # Original plane: keep all elevations positive because the current geometry
    # initialization clips negative DEM elevations to zero.
    slope = np.tan(np.deg2rad(args.slope_deg))
    z_ref = 10.0
    B0 = z_ref - slope * Xc

    # Rectangular excavation, finite in both x and y.
    x1, x2 = 10.0, 20.0
    y1, y2 = -5.0, 5.0
    depth = args.depth

    mask = (Xc >= x1) & (Xc < x2) & (Yc >= y1) & (Yc < y2)

    B = B0.copy()
    B[mask] -= depth

    H = np.zeros_like(B)
    H[mask] = depth

    U = np.zeros_like(B)
    V = np.zeros_like(B)
    eta = B + H

    # Exact intended relation in the wet region: eta = original plane B0.
    err = np.max(np.abs(eta[mask] - B0[mask])) if np.any(mask) else np.nan

    # Keep the test geometrically well resolved and the nominal step exactly at
    # computational interfaces whenever the default dimensions are used.
    print(f"nx_cells = {nx_cells}")
    print(f"ny_cells = {ny_cells}")
    print(f"dx = dy = {dx:.12g} m")
    print(f"slope angle = {args.slope_deg:.6g} deg")
    print(f"excavation depth = {depth:.6g} m")
    print(f"excavation: x=[{x1},{x2}), y=[{y1},{y2})")
    print(f"max |eta-B0| in wet cells = {err:.3e} m")
    print(f"initial water volume = {np.sum(H) * dx * dx:.12g} m^3")

    if np.min(B) < 0.0:
        raise ValueError(
            "Topography became negative. Increase z_ref or reduce slope/depth; "
            "the current code clips negative DEM elevations to zero."
        )

    outdir = Path(__file__).resolve().parent
    topo_file = outdir / "topography_dem.asc"
    write_esri_cell_centered(
        topo_file, B, x_min, y_min, dx, slope, z_ref, x1, x2, y1, y2, depth
    )

    # Restart format follows EXAMPLE_2D / EXAMPLE_BUMP conventions. With
    # N_SOLID=0, N_ADD_GAS=0, LIQUID_FLAG=T, we write 6 columns:
    # x, y, rho*h, rho*h*u, rho*h*v, rho*h*(cp*T + kinetic energy).
    init_file = outdir / "inclined_excavation_2D_0000.q_2d"
    with init_file.open("w") as f:
        for j in range(ny_cells):
            q0 = np.zeros((6, nx_cells))
            q0[0, :] = Xc[j, :]
            q0[1, :] = Yc[j, :]
            q0[2, :] = rho_l * H[j, :]
            q0[3, :] = rho_l * H[j, :] * U[j, :]
            q0[4, :] = rho_l * H[j, :] * V[j, :]
            q0[5, :] = rho_l * H[j, :] * (
                sp_heat_l * temperature + 0.5 * (U[j, :]**2 + V[j, :]**2)
            )
            np.savetxt(f, q0.T, fmt="%19.12e")
            f.write(" \n")

    template = (outdir / "IMEX_SfloW2D.template").read_text()
    inp = template.replace("runname", "inclinedExcavation2D")
    inp = inp.replace("restartfile", init_file.name)
    inp = inp.replace("x_min", str(x_min))
    inp = inp.replace("y_min", str(y_min))
    inp = inp.replace("nx_cells", str(nx_cells))
    inp = inp.replace("ny_cells", str(ny_cells))
    inp = inp.replace("dx", str(dx))
    (outdir / "IMEX_SfloW2D.inp").write_text(inp)

    # Save compact diagnostics for regression tests.
    np.savez(
        outdir / "initial_condition_reference.npz",
        x=x_cent,
        y=y_cent,
        B0=B0,
        B=B,
        H=H,
        eta=eta,
        mask=mask,
        dx=dx,
        x1=x1,
        x2=x2,
        y1=y1,
        y2=y2,
        slope_deg=args.slope_deg,
        depth=depth,
    )

    if args.plot:
        import matplotlib.pyplot as plt

        # Plan-view thickness.
        fig, ax = plt.subplots(figsize=(8, 4.6))
        m = ax.pcolormesh(Xc, Yc, H, shading="nearest")
        fig.colorbar(m, ax=ax, label="h [m]")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title("Initial water thickness")
        fig.tight_layout()
        fig.savefig(outdir / "initial_thickness.png", dpi=180)
        plt.close(fig)

        # Centerline longitudinal section.
        j0 = int(np.argmin(np.abs(y_cent)))
        fig, ax = plt.subplots(figsize=(8, 4.6))
        ax.plot(x_cent, B0[j0, :], label="original plane B0")
        ax.plot(x_cent, B[j0, :], label="excavated bed B")
        ax.plot(x_cent, eta[j0, :], label="free surface eta")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("elevation [m]")
        ax.set_title(f"Centerline section (y={y_cent[j0]:.3f} m)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "initial_centerline.png", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    main()
