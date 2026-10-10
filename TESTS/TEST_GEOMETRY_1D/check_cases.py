"""Independent refreshed-geometry oracle, including degenerate and short lines."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np

HERE = Path(__file__).resolve().parent
CONTRACT = json.loads((HERE / "contract.json").read_text())
TOL = CONTRACT["roundoff_epsilon_multiplier"] * np.finfo(float).eps


def expected(nx, ny, kind, slope):
    """Differentiate independently fitted polynomials; retain the specified 2D kernel."""
    dx, dy = CONTRACT["spacing"]
    x, y = np.meshgrid(np.arange(nx+1)*dx, np.arange(ny+1)*dy)
    if nx == 1:
        x[:] = 0
    if ny == 1:
        y[:] = 0
    vertices = np.full(x.shape, 10.0)
    if kind >= 1:
        vertices += .15*x-.11*y
    if kind >= 2:
        vertices += .02*x*x-.01*y*y+.013*x*y
    if kind == 3:
        vertices += .07*np.sin(.7*x+.3*y)
    fx = .5*(vertices[:-1]+vertices[1:])
    fy = .5*(vertices[:, :-1]+vertices[:, 1:])
    center = .25*(vertices[:-1, :-1]+vertices[:-1, 1:]+vertices[1:, :-1]+vertices[1:, 1:])
    derivatives = np.zeros((5, ny, nx))
    if min(nx, ny) == 1:
        line = center[0] if ny == 1 else center[:, 0]
        spacing = dx if ny == 1 else dy
        count = len(line)
        first, second = np.zeros(count), np.zeros(count)
        if count > 1:
            coordinates = (np.arange(count)-(count-1)/2)*spacing
            if count < 5:
                coefficients = np.polynomial.polynomial.polyfit(coordinates, line, min(2, count-1))
                first[:] = coefficients[1]
                if count >= 3:
                    first += 2*coefficients[2]*coordinates
                    second[:] = 2*coefficients[2]
            else:
                for i in range(2, count-2):
                    coefficients = np.polynomial.polynomial.polyfit(np.arange(-2, 3)*spacing, line[i-2:i+3], 2)
                    first[i] = coefficients[1]; second[i] = 2*coefficients[2]
                first[:2] = first[2]; first[-2:] = first[-3]
                second[:2] = second[2]; second[-2:] = second[-3]
        if ny == 1:
            derivatives[0, 0] = first; derivatives[2, 0] = second
        else:
            derivatives[1, :, 0] = first; derivatives[3, :, 0] = second
    else:
        c1, c2 = np.array([-2, -1, 0, 1, 2]), np.array([2, -1, -2, -1, 2])
        for k in range(2, ny-2):
            for j in range(2, nx-2):
                derivatives[:, k, j] = [c1@center[k, j-2:j+3]/(10*dx),
                    c1@center[k-2:k+3, j]/(10*dy), c2@center[k, j-2:j+3]/(7*dx*dx),
                    c2@center[k-2:k+3, j]/(7*dy*dy),
                    np.sum(np.outer(c1, c1)*center[k-2:k+3, j-2:j+3])/(100*dx*dy)]
        derivatives[:, :, :2] = derivatives[:, :, 2:3]
        derivatives[:, :, -2:] = derivatives[:, :, -3:-2]
        derivatives[:, :2] = derivatives[:, 2:3]
        derivatives[:, -2:] = derivatives[:, -3:-2]
    g = 1/(1+derivatives[0]**2+derivatives[1]**2) if slope else np.ones_like(center)
    gx = np.c_[g[:, 0], .5*(g[:, :-1]+g[:, 1:]), g[:, -1]]
    gy = np.r_[g[:1], .5*(g[:-1]+g[1:]), g[-1:]]
    return [vertices, fx, fy, center, *derivatives, g, gx, gy]


def run(executable, profile, baseline_only=False):
    """Require raw byte equality and independent field agreement for every frozen case."""
    records = []
    grids = [(n, 1) for n in CONTRACT["lengths"]] + [(1, n) for n in CONTRACT["lengths"]]
    grids += [tuple(CONTRACT["control_2d"])]
    grids = list(dict.fromkeys(grids))  # The single-cell grid is common to both axes.
    if baseline_only:
        grids = [(nx, ny) for nx, ny in grids if nx >= 5 and (ny == 1 or ny >= 5)]
    for nx, ny in grids:
        for bed in CONTRACT["beds"]:
            for slope in CONTRACT["slope_flags"]:
                label = f"grid-{nx}-{ny}-bed-{bed}-S{int(slope)}"
                snapshots = []
                for threads in (1, 4):
                    work = Path(label)/f"threads-{threads}"; work.mkdir(parents=True)
                    env = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE"}
                    result = subprocess.run([str(executable), str(threads), str(nx), str(ny), str(bed), str(int(slope))],
                        cwd=work, env=env, capture_output=True, text=True, timeout=30)
                    (work/"solver.log").write_text(result.stdout+result.stderr)
                    if result.returncode:
                        raise AssertionError(label+": "+result.stdout+result.stderr)
                    snapshots.append((work/"geometry.bin").read_bytes())
                if snapshots[0] != snapshots[1]:
                    raise AssertionError(label+": actual 1/4-thread mismatch")
                raw = np.frombuffer(snapshots[0], dtype=np.float64); cursor = 0; error = 0
                for target in expected(nx, ny, bed, slope):
                    observed = raw[cursor:cursor+target.size].reshape(target.shape)
                    error = max(error, float(np.max(np.abs(observed-target)/np.maximum(1, np.abs(target)))))
                    cursor += target.size
                if cursor != raw.size or not np.all(np.isfinite(raw)) or error > TOL:
                    raise AssertionError(f"{label}: independent geometry error {error}, limit {TOL}")
                records.append({"case": label, "error": error, "actual_threads": [1, 4],
                    "sha256": hashlib.sha256(snapshots[0]).hexdigest()})
    if len(records) != (32 if baseline_only else 112):
        raise AssertionError("incomplete geometry inventory")
    Path("evidence.json").write_text(json.dumps({"profile": profile, "contract": CONTRACT,
        "baseline_subset": baseline_only, "cases": records}, indent=2)+"\n")
    print(f"PASS: {profile} geometry, {len(records)} cases, actual one/four teams")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("profile")
    parser.add_argument("--baseline-only", action="store_true")
    args = parser.parse_args()
    run(args.executable.resolve(), args.profile, args.baseline_only)
