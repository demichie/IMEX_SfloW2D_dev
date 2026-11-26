#!/usr/bin/env python3
"""Plot ESRI ASCII solver rasters over the exported computational terrain."""
import argparse
import csv
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LightSource


def read_asc(path):
    """Return the raster and its cell-centre geographic extent."""
    path = Path(path)
    with path.open(encoding="utf-8") as stream:
        fields = {}
        for _ in range(6):
            key, value = stream.readline().split()
            fields[key.lower()] = float(value)
    cols, rows = int(fields["ncols"]), int(fields["nrows"])
    cell = fields["cellsize"]
    if cols < 1 or rows < 1 or cell <= 0:
        raise ValueError(f"Invalid raster dimensions: {path}")
    x = fields.get("xllcorner")
    y = fields.get("yllcorner")
    if x is None:
        x = fields["xllcenter"] - cell / 2
    if y is None:
        y = fields["yllcenter"] - cell / 2
    arr = np.loadtxt(path, skiprows=6, ndmin=2)
    if arr.shape != (rows, cols):
        raise ValueError(f"Unexpected raster shape in {path}: {arr.shape}")
    nodata = fields.get("nodata_value")
    if nodata is not None:
        arr[arr == nodata] = np.nan
    return arr, (x, x + cols * cell, y, y + rows * cell)


def manifest_times(run_dir):
    path = run_dir / "metadata" / "output_manifest.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ["snapshot_index", "time_seconds"]:
            raise ValueError(f"Invalid snapshot manifest: {path}")
        times = {}
        for row in reader:
            index = int(row["snapshot_index"])
            if index in times:
                raise ValueError(f"Duplicate snapshot index: {index}")
            times[index] = float(row["time_seconds"])
    return times


def plot_raster(source, terrain, extent, target, label):
    raster, raster_extent = read_asc(source)
    if raster.shape != terrain.shape or not np.allclose(raster_extent, extent):
        raise ValueError(f"Raster grid differs from terrain: {source}")
    if not np.any(np.isfinite(raster)):
        return False
    cell = (extent[1] - extent[0]) / raster.shape[1]
    hillshade = LightSource(azdeg=315, altdeg=45).hillshade(
        np.nan_to_num(terrain, nan=0.0), dx=cell, dy=cell)
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.imshow(hillshade, extent=extent, origin="upper", cmap="gray")
    layer = ax.imshow(np.ma.masked_invalid(raster), extent=extent,
                      origin="upper", cmap="viridis", alpha=0.7)
    fig.colorbar(layer, ax=ax, label=label)
    ax.set(xlabel="X (m)", ylabel="Y (m)")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(target, dpi=160)
    plt.close(fig)
    return True


def run(run_dir, kind="all"):
    run_dir = Path(run_dir).resolve()
    if not run_dir.is_dir():
        raise FileNotFoundError(f"Run output directory not found: {run_dir}")
    terrain, extent = read_asc(run_dir / "terrain" / "dem_esri.asc")
    written = []
    if kind in ("all", "snapshots"):
        times = manifest_times(run_dir)
        for directory in sorted((run_dir / "snapshots").glob("*_t_*")):
            if not directory.is_dir():
                continue
            index = int(directory.name.split("_t_", 1)[0])
            if index not in times:
                raise ValueError(f"Unlisted snapshot directory: {directory}")
            for source in sorted(directory.glob("*.asc")):
                label = f"{source.stem} at t={times[index]:.6g} s"
                target = source.with_suffix(".png")
                if plot_raster(source, terrain, extent, target, label):
                    written.append(target)
    if kind in ("all", "diagnostics"):
        for source in sorted((run_dir / "diagnostics").glob("*.asc")):
            target = source.with_suffix(".png")
            if plot_raster(source, terrain, extent, target, source.stem):
                written.append(target)
    if kind in ("all", "vulnerability"):
        for source in sorted((run_dir / "vulnerability").glob("*_VT_*.asc")):
            target = source.with_suffix(".png")
            if plot_raster(source, terrain, extent, target, source.stem):
                written.append(target)
    return written


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path, help="OUTPUT/<RUN_NAME> directory")
    parser.add_argument("--kind", choices=("all", "snapshots", "diagnostics",
                                           "vulnerability"), default="all")
    args = parser.parse_args()
    for path in run(args.run_dir, args.kind):
        print(path)


if __name__ == "__main__":
    main()
