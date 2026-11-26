#!/usr/bin/env python3
"""Etna two-layer integration smoke test for output and restart files."""
import argparse
import csv
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from prepare import EXAMPLE, FIELD_WIDTH, N_VARS, prepare


def change_parameter(text, name, value):
    pattern = re.compile(r"(?im)^(\s*" + re.escape(name) + r"\s*=\s*)[^,\n]+(,?)")
    updated, count = pattern.subn(lambda match: match.group(1) + value + ",", text)
    if count != 1:
        raise ValueError(f"Expected exactly one {name}, found {count}")
    return updated


def check_snapshot(path, expect_dry=False):
    if not path.is_file():
        raise FileNotFoundError(f"Missing conservative snapshot: {path}")
    cells = 0
    total_mass = 0.0
    with path.open(encoding="ascii") as stream:
        for line in stream:
            if not line.strip():
                continue
            record = line.rstrip("\r\n")
            if len(record) != FIELD_WIDTH * (N_VARS + 2):
                raise ValueError(f"Invalid fixed-width snapshot record in {path}")
            values = [float(record[i:i + FIELD_WIDTH].replace("D", "E"))
                      for i in range(0, len(record), FIELD_WIDTH)]
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"Non-finite conservative field in {path}")
            if values[2] < -1e-12:
                raise ValueError(f"Negative layer mass in {path}")
            if expect_dry and any(value != 0.0 for value in values[2:]):
                raise ValueError(f"Upper layer not initially dry: {path}")
            total_mass += values[2]
            cells += 1
    if cells == 0:
        raise ValueError(f"Empty conservative snapshot: {path}")
    return cells, total_mass


def check_run(root, run_name):
    out = root / "OUTPUT" / run_name
    manifest = out / "metadata" / "output_manifest.csv"
    with manifest.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) < 2 or [int(row["snapshot_index"]) for row in rows[:2]] != [0, 1]:
        raise ValueError("Expected initial and advanced output snapshots")
    times = [float(row["time_seconds"]) for row in rows]
    if not all(math.isfinite(t) for t in times) or any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError("Snapshot times are not finite and strictly increasing")

    for index in (0, 1):
        matching = list((out / "snapshots").glob(f"{index:04d}_t_*"))
        if len(matching) != 1:
            raise ValueError(f"Expected one snapshot directory for index {index}")
        snap = matching[0]
        l1 = snap / f"{run_name}_L1_{index:04d}.q_2d"
        l2 = snap / f"{run_name}_L2_{index:04d}.q_2d"
        cells1, mass1 = check_snapshot(l1)
        cells2, mass2 = check_snapshot(l2, expect_dry=(index == 0))
        if cells1 != cells2 or mass1 <= 0:
            raise ValueError("Mismatched layer grids or nonpositive lower-layer mass")
        for suffix in (f"_L2_{index:04d}.asc", f"_T_L2_{index:04d}.asc",
                       f"_total_depth_{index:04d}.asc", f"_eta_L2_{index:04d}.asc"):
            if not (snap / (run_name + suffix)).is_file():
                raise FileNotFoundError(f"Missing two-layer raster: {run_name + suffix}")
        print(f"Snapshot {index}: time={times[index]:.8g}, L1 mass={mass1:.8g}, L2 mass={mass2:.8g}")

    if not (out / "restart" / "restart.bin").is_file():
        raise FileNotFoundError("Missing canonical binary restart")
    if not (out / "fields" / (run_name + ".nc")).is_file():
        raise FileNotFoundError("Missing canonical NetCDF output")
    budgets = sorted((out / "fields").glob(run_name + "_budget_*.csv"))
    if len(budgets) < 2:
        raise ValueError("Expected per-snapshot two-layer budgets")
    for budget in budgets:
        with budget.open(newline="", encoding="utf-8") as stream:
            reader = csv.reader(stream)
            header = next(reader)
            values = [float(v) for v in next(reader)]
        if "mass_L1" not in ",".join(header) or not all(map(math.isfinite, values)):
            raise ValueError(f"Invalid per-layer budget: {budget}")
        if len(values) < 4 or abs(values[3] - values[1] - values[2]) > (
            1e-9 * max(1.0, abs(values[3]))
        ):
            raise ValueError(f"Inconsistent L1+L2 mass budget: {budget}")
    print("PASS: two-layer initialization, advanced output, CSV budgets, NetCDF and restart files")


def run_test(solver, root, end_time, timeout):
    root.mkdir(parents=True, exist_ok=True)
    text = (EXAMPLE / "IMEX_SfloW2D.inp").read_text(encoding="utf-8")
    (root / "IMEX_SfloW2D.inp").write_text(text, encoding="utf-8")
    prepare(root)
    text = change_parameter(text, "T_END", f"{end_time:.16g}")
    text = change_parameter(text, "DT_OUTPUT", f"{end_time:.16g}")
    (root / "IMEX_SfloW2D.inp").write_text(text, encoding="utf-8")
    env = dict(os.environ, OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE")
    process = subprocess.run([str(solver)], cwd=root, env=env, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             timeout=timeout, check=False)
    (root / "runner.log").write_text(process.stdout, encoding="utf-8")
    if process.returncode != 0:
        print(process.stdout[-16000:], file=sys.stderr)
        raise RuntimeError(f"Two-layer solver exited with status {process.returncode}")
    check_run(root, "exampleEtna_two_layer")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", required=True, type=Path, help="Compiled solver executable")
    parser.add_argument("--workdir", type=Path, help="New empty directory; retained for inspection")
    parser.add_argument("--end-time", type=float, default=0.001)
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args()
    if not 0 < args.end_time <= 0.1:
        parser.error("Smoke-test end time must be between zero and 0.1 s")
    solver = args.solver.resolve()
    if not solver.is_file():
        parser.error(f"Solver binary not found: {solver}")
    if args.workdir:
        root = args.workdir.resolve()
        if root.exists() and any(root.iterdir()):
            parser.error("Work directory must be absent or empty")
        run_test(solver, root, args.end_time, args.timeout)
    else:
        with tempfile.TemporaryDirectory(prefix="imex-two-layer-etna-") as directory:
            run_test(solver, Path(directory), args.end_time, args.timeout)


if __name__ == "__main__":
    main()
