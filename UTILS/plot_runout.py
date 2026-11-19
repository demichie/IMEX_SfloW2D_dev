#!/usr/bin/env python3
"""Plot solver runout and area history, reporting maximum-runout milestones."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def analyze(path):
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        next(reader, None)
        rows = [row for row in reader if any(part.strip() for part in row)]
    if not rows:
        raise ValueError(f"No runout samples in {path}")
    data = np.atleast_2d(np.asarray([[float(v) for v in row] for row in rows]))
    if data.shape[1] < 3 or not np.all(np.isfinite(data[:, :3])):
        raise ValueError(f"Invalid runout CSV: {path}")
    time, distance, area = data[:, 0], data[:, 1], data[:, 2]
    if np.any(np.diff(time) <= 0):
        raise ValueError("Runout times must be strictly increasing")
    fig, ax = plt.subplots(figsize=(9, 5))
    twin = ax.twinx()
    ax.plot(time, distance, "C0-", label="Runout")
    twin.plot(time, area, "C1-", label="Area")
    ax.set(xlabel="Time (s)", ylabel="Runout (m)")
    twin.set_ylabel("Area (m²)")
    fig.tight_layout()
    output = path.with_suffix(".png")
    fig.savefig(output)
    plt.close(fig)
    maximum = float(np.max(distance))
    milestones = {}
    if maximum > 0:
        for portion in (0.95, 0.99):
            indices = np.flatnonzero(distance >= portion * maximum)
            milestones[portion] = (float(time[indices[0]]), float(distance[indices[0]]))
    return output, float(time[np.argmax(distance)]), maximum, milestones


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runout_csv", type=Path)
    args = parser.parse_args()
    output, when, maximum, milestones = analyze(args.runout_csv)
    print(f"{output}\nMaximum runout: {maximum:g} m at {when:g} s")
    for fraction, (time, distance) in milestones.items():
        print(f"First {fraction:.0%} of maximum: {distance:g} m at {time:g} s")


if __name__ == "__main__":
    main()
