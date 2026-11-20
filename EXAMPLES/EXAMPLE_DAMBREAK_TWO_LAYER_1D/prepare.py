#!/usr/bin/env python3
"""Generate a genuinely one-cell-wide flat-bed two-layer dam-break test.

Only input rasters are generated; the Fortran solver is never launched.
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
NX = 200
NY = 1
DX = 1.0
DAM_CELLS = 60
DAM_HEIGHT = 2.0


def write_asc(path, ncols, nrows, rows):
    header = (
        f"ncols {ncols}\n"
        f"nrows {nrows}\n"
        "xllcorner 0\n"
        "yllcorner 0\n"
        f"cellsize {DX:g}\n"
        "NODATA_value -9999\n"
    )
    with path.open("w", encoding="ascii") as stream:
        stream.write(header)
        for row in rows:
            stream.write(" ".join(f"{value:.8f}" for value in row) + "\n")


def main():
    dem = HERE / "DEM"
    dem.mkdir(parents=True, exist_ok=True)
    write_asc(dem / "topography.asc", NX + 1, NY + 1,
              [[0.0] * (NX + 1) for _ in range(NY + 1)])
    write_asc(dem / "l1_thickness.asc", NX, NY,
              [[DAM_HEIGHT if j < DAM_CELLS else 0.0 for j in range(NX)]])
    print(f"Prepared 1D dam-break: {NX} x {NY} cells, dx={DX} m, "
          f"initial L1 volume={DAM_CELLS * DAM_HEIGHT * DX**2:g} m3.")
    print("L2 is initially dry; the solver was not executed.")


if __name__ == "__main__":
    main()
