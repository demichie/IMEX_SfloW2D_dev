#!/usr/bin/env python3
"""TL13a extended two-layer Etna checks: evolution, restart and OpenMP.

This is functional regression evidence, not validation of coupled eigenvalues
or an endorsement of two-layer physics for production use.
"""

import argparse
import csv
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from prepare import EXAMPLE, FIELD_WIDTH, N_VARS, NX, NY, prepare
from run_validation import change_parameter, check_run

RUN = "exampleEtna_two_layer"
BUDGET_COLUMNS = (
    "time", "mass_L1", "mass_L2", "mass_total", "volume_L1",
    "volume_L2", "volume_total", "solid_L1", "solid_L2",
    "solid_substrate", "solid_system",
)
REL_CONSERVATION = 1.0e-6
REL_COMPARISON = 5.0e-10
ABS_COMPARISON = 1.0e-8
REL_SNAPSHOT_BUDGET = 2.0e-8


def set_input(root, end_time, output_interval, restart_binary=False):
    """Keep the canonical two-layer example; alter only the run horizon."""
    text = (EXAMPLE / "IMEX_SfloW2D.inp").read_text(encoding="utf-8")
    for key, value in (
        ("T_END", f"{end_time:.16g}"),
        ("DT_OUTPUT", f"{output_interval:.16g}"),
        ("OUTPUT_RUNOUT_FLAG", "F"),
    ):
        text = change_parameter(text, key, value)
    if restart_binary:
        text = change_parameter(text, "N_RESTART_FILES", "1")
        original = re.compile(r"(?im)^[ \t]*RESTART_FILES[ \t]*=.*$")
        checkpoint = f"'OUTPUT/{RUN}/restart/restart.bin',"
        text, count = original.subn(" RESTART_FILES=" + checkpoint, text)
        if count != 1:
            raise ValueError("Expected exactly one RESTART_FILES assignment")
    (root / "IMEX_SfloW2D.inp").write_text(text, encoding="utf-8")


def execute(solver, root, threads, timeout, label):
    env = dict(os.environ, IMEX_TWO_LAYER_VALIDATION="1",
               OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1")
    try:
        result = subprocess.run([str(solver)], cwd=root, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                check=False, timeout=timeout)
    except subprocess.TimeoutExpired as error:
        raise AssertionError(f"{label}: timed out after {timeout} s") from error
    (root / (label + ".log")).write_text(result.stdout, encoding="utf-8")
    if result.returncode:
        raise AssertionError(
            f"{label}: solver exit {result.returncode}\n" + result.stdout[-12000:]
        )
    match = re.search(r"Parallel run:\s*number of threads used\s+(\d+)", result.stdout)
    if match is None or int(match.group(1)) != threads:
        raise AssertionError(f"{label}: OpenMP team size {threads} not verified")
    if "TL13a experimental two-layer validation" not in result.stdout:
        raise AssertionError(f"{label}: validation gate was not acknowledged")


def budgets(root):
    run_root = root / "OUTPUT" / RUN
    manifest = run_root / "metadata" / "output_manifest.csv"
    with manifest.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ["snapshot_index", "time_seconds"]:
            raise AssertionError("Unexpected output manifest schema")
        entries = list(reader)
    if len(entries) < 3:
        raise AssertionError("Extended trajectory must contain at least three output times")
    times = [float(entry["time_seconds"]) for entry in entries]
    indices = [int(entry["snapshot_index"]) for entry in entries]
    if indices != list(range(len(entries))) or any(
        not math.isfinite(t) for t in times
    ) or any(b <= a for a, b in zip(times, times[1:])):
        raise AssertionError("Invalid output time/index sequence")
    result = []
    for index, time in enumerate(times):
        path = run_root / "fields" / f"{RUN}_budget_{index:04d}.csv"
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if tuple(reader.fieldnames or ()) != BUDGET_COLUMNS:
                raise AssertionError(f"Unexpected budget schema in {path}")
            rows = list(reader)
        if len(rows) != 1:
            raise AssertionError(f"Expected one budget record in {path}")
        record = {key: float(rows[0][key]) for key in BUDGET_COLUMNS}
        if not all(math.isfinite(value) for value in record.values()):
            raise AssertionError(f"Nonfinite two-layer budget in {path}")
        if abs(record["time"] - time) > 1e-12 * max(1.0, abs(time)):
            raise AssertionError("Manifest/budget times disagree")
        for total, parts in (
            ("mass_total", ("mass_L1", "mass_L2")),
            ("volume_total", ("volume_L1", "volume_L2")),
            ("solid_system", ("solid_L1", "solid_L2", "solid_substrate")),
        ):
            expected = sum(record[key] for key in parts)
            if abs(record[total] - expected) > 1e-10 * max(1.0, abs(expected)):
                raise AssertionError(f"{total} inconsistent at index {index}")
        for key in ("mass_L1", "mass_L2", "solid_L1", "solid_L2",
                    "solid_substrate", "volume_L1", "volume_L2"):
            if record[key] < -1e-10:
                raise AssertionError(f"Negative {key} at index {index}")
        result.append(record)
    return result


def snapshot(root, index, layer):
    matches = list((root / "OUTPUT" / RUN / "snapshots").glob(f"{index:04d}_t_*"))
    if len(matches) != 1:
        raise AssertionError(f"Snapshot {index} not unique in {root}")
    path = matches[0] / f"{RUN}_L{layer}_{index:04d}.q_2d"
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def cells(path):
    with path.open(encoding="ascii") as stream:
        for line in stream:
            if not line.strip():
                continue
            row = line.rstrip("\r\n")
            if len(row) != FIELD_WIDTH * (N_VARS + 2):
                raise AssertionError(f"Invalid conservative snapshot record in {path}")
            values = [float(row[i:i + FIELD_WIDTH].replace("D", "E"))
                      for i in range(0, len(row), FIELD_WIDTH)]
            if not all(math.isfinite(value) for value in values):
                raise AssertionError(f"Nonfinite conservative snapshot in {path}")
            if values[2] < -1e-10 or values[6] < -1e-10:
                raise AssertionError(f"Negative total/solid inventory in {path}")
            yield values


def compare_snapshots(reference, candidate, index, tolerance=REL_COMPARISON):
    for layer in (1, 2):
        a = cells(snapshot(reference, index, layer))
        b = cells(snapshot(candidate, index, layer))
        count = 0
        from itertools import zip_longest
        for left, right in zip_longest(a, b):
            if left is None or right is None:
                raise AssertionError("Different snapshot cell counts")
            for value, other in zip(left, right):
                if abs(value - other) > ABS_COMPARISON + tolerance * abs(value):
                    raise AssertionError(
                        f"Snapshot index {index}, L{layer}, cell {count}: "
                        f"{value:.16g} != {other:.16g}"
                    )
            count += 1
        if count == 0:
            raise AssertionError("Empty comparison snapshot")


def check_invariants(root, require_growth=True):
    records = budgets(root)
    for key in ("mass_total", "solid_system"):
        initial = records[0][key]
        for n, record in enumerate(records):
            if abs(record[key] - initial) > REL_CONSERVATION * max(1.0, abs(initial)):
                raise AssertionError(
                    f"Conservation failure {key} at output {n}: "
                    f"{record[key]:.16g} vs {initial:.16g}"
                )
    if records[0]["mass_L2"] != 0.0 or records[0]["solid_L2"] != 0.0:
        raise AssertionError("L2 must be exactly dry at initial output")
    max_mass_l2 = max(record["mass_L2"] for record in records[1:])
    if require_growth and max_mass_l2 <= 0.0:
        raise AssertionError("No measurable L2 growth during extended elutriation run")
    input_text = (root / "IMEX_SfloW2D.inp").read_text(encoding="utf-8")
    match = re.search(r"(?im)^[ \t]*CELL_SIZE[ \t]*=[ \t]*([0-9.+EeDd-]+)", input_text)
    if match is None:
        raise AssertionError("Missing canonical cell size")
    cell_size = float(match.group(1).replace("D", "E").replace("d", "e"))
    if not math.isfinite(cell_size) or cell_size <= 0.0:
        raise AssertionError("Invalid canonical cell size")
    cell_area = cell_size * cell_size
    for n, record in enumerate(records):
        for layer in (1, 2):
            count = 0
            mass = 0.0
            solids = 0.0
            for state in cells(snapshot(root, n, layer)):
                mass += state[2]
                solids += state[6]
                count += 1
            if count != NX * NY:
                raise AssertionError(f"Unexpected snapshot grid size {count}")
            for key, actual in ((f"mass_L{layer}", mass * cell_area),
                                (f"solid_L{layer}", solids * cell_area)):
                expected = record[key]
                if abs(actual - expected) > REL_SNAPSHOT_BUDGET * max(1.0, abs(expected)):
                    raise AssertionError(
                        f"Conservative snapshot/{key} mismatch at index {n}: "
                        f"{actual:.16g} vs {expected:.16g}"
                    )
    print(f"PASS {root.name}: {len(records)} snapshots, "
          f"max L2 mass={max_mass_l2:.12g}, conserved mass and solids")
    return records


def compare_budgets(reference, candidate):
    left, right = budgets(reference), budgets(candidate)
    if len(left) != len(right):
        raise AssertionError("Different output counts")
    for n, (a, b) in enumerate(zip(left, right)):
        for key in BUDGET_COLUMNS:
            if abs(a[key] - b[key]) > ABS_COMPARISON + REL_COMPARISON * abs(a[key]):
                raise AssertionError(f"Budget mismatch index {n} field {key}")


def run_validation(solver, root, horizon, timeout):
    midpoint = horizon / 2.0
    template = root / "initial"
    template.mkdir(parents=True)
    prepare(template, solver)
    runs = {name: root / name for name in ("reference", "resumed", "threads2", "threads4")}
    for directory in runs.values():
        shutil.copytree(template, directory)
        set_input(directory, horizon, midpoint)
    reference = runs["reference"]
    execute(solver, reference, 1, timeout, "continuous")
    check_run(reference, RUN)
    reference_budgets = check_invariants(reference)
    if len(reference_budgets) != 3:
        raise AssertionError("Expected exactly three scheduled snapshots")

    resumed = runs["resumed"]
    set_input(resumed, midpoint, midpoint)
    execute(solver, resumed, 1, timeout, "first_segment")
    first_checkpoint = resumed / "OUTPUT" / RUN / "restart" / "restart.bin"
    if not first_checkpoint.is_file():
        raise AssertionError("Missing intermediate binary restart")
    set_input(resumed, horizon, midpoint, restart_binary=True)
    execute(solver, resumed, 1, timeout, "second_segment")
    check_run(resumed, RUN)
    check_invariants(resumed)
    compare_budgets(reference, resumed)
    compare_snapshots(reference, resumed, 2)
    if (reference / "OUTPUT" / RUN / "restart" / "restart.bin").read_bytes() != (
        resumed / "OUTPUT" / RUN / "restart" / "restart.bin"
    ).read_bytes():
        raise AssertionError("Final binary checkpoints differ after restart continuation")
    print("PASS: continuous/segmented two-layer restart and final checkpoint")

    for threads in (2, 4):
        candidate = runs[f"threads{threads}"]
        execute(solver, candidate, threads, timeout, f"threads{threads}")
        check_run(candidate, RUN)
        check_invariants(candidate)
        compare_budgets(reference, candidate)
        for index in (1, 2):
            compare_snapshots(reference, candidate, index)
        print(f"PASS: two-layer OpenMP {threads} threads vs 1 thread")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True,
                        help="New or empty evidence directory")
    parser.add_argument("--horizon", type=float, default=0.01,
                        help="Final simulation time in seconds (default 0.01)")
    parser.add_argument("--timeout", type=int, default=900,
                        help="Per-run timeout in seconds")
    args = parser.parse_args()
    if not 0.0 < args.horizon <= 0.1:
        parser.error("The validation horizon must satisfy 0 < t <= 0.1 s")
    root = args.workdir.resolve()
    if root.exists() and any(root.iterdir()):
        parser.error("Evidence directory must be new or empty")
    root.mkdir(parents=True, exist_ok=True)
    try:
        run_validation(args.solver.resolve(), root, args.horizon, args.timeout)
    except Exception as error:
        print(f"FAIL: TL13a extended validation: {error}", file=sys.stderr)
        return 1
    print("PASS: extended TL13a functional validation (not scientific approval)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
