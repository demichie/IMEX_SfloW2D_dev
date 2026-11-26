#!/usr/bin/env python3
"""Plot each time series from a solver probe CSV, optionally with its derivative."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(path, rate=False):
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        coordinates = next(reader, None)
        names = next(reader, None)
        units = next(reader, None)
        if not coordinates or not names or not units or len(names) != len(units):
            raise ValueError(f"Invalid probe CSV headers: {path}")
        rows = [row for row in reader if any(part.strip() for part in row)]
    if not rows:
        raise ValueError(f"No probe samples in {path}")
    data = np.atleast_2d(np.asarray([[float(v) for v in row] for row in rows]))
    if data.shape[1] != len(names):
        raise ValueError(f"Unexpected number of probe columns: {path}")
    time = data[:, 0]
    if not np.all(np.isfinite(time)) or np.any(np.diff(time) <= 0):
        raise ValueError("Probe times must be finite and strictly increasing")
    outputs = []
    for col in range(1, data.shape[1]):
        name = names[col].strip()
        unit = units[col].strip()
        safe_name = "".join(ch if ch.isalnum() or ch in "_-" else "_" for ch in name)
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(time, data[:, col])
        ax.set(xlabel="Time (s)", ylabel=f"{name} {unit}")
        target = path.with_name(f"{path.stem}_{safe_name}.pdf")
        fig.tight_layout()
        fig.savefig(target)
        outputs.append(target)
        if rate:
            if len(time) < 2:
                raise ValueError("Two or more samples required for a rate plot")
            slope = np.diff(data[:, col]) / np.diff(time)
            second = ax.twinx()
            second.plot(0.5 * (time[1:] + time[:-1]), slope, "C1--")
            second.set_ylabel(f"Rate of {name} ({unit}/s)")
            rate_target = path.with_name(f"{path.stem}_{safe_name}_rate.pdf")
            fig.tight_layout()
            fig.savefig(rate_target)
            outputs.append(rate_target)
        plt.close(fig)
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probe_csv", type=Path)
    parser.add_argument("--rate", action="store_true")
    args = parser.parse_args()
    for output in plot(args.probe_csv, args.rate):
        print(output)


if __name__ == "__main__":
    main()
