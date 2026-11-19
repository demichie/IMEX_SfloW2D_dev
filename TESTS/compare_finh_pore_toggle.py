#!/usr/bin/env python3
"""Compare the unchanged TEST_OFF_FINH input with a no-pore-pressure copy."""

import argparse
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time

from compare_finh_timestep import extract_timesteps

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "TESTS" / "TEST_DAMBREAK" / "TEST_OFF_FINH"
FLAG = re.compile(r"^(\s*PORE_PRESSURE_FLAG\s*=\s*)T(\s*,?\s*)$", re.I | re.M)
GROUP = re.compile(r"^\s*&PORE_PRESSURE_PARAMETERS\b[^\n]*\n(?:[^\n]*\n)*?\s*/\s*(?:\n|$)", re.I | re.M)


def without_pore_pressure(text):
    """Only alter the temporary input; never modify the tracked fixture."""
    text, flags = FLAG.subn(lambda match: match.group(1) + "F" + match.group(2), text)
    text, groups = GROUP.subn("", text)
    if flags != 1 or groups != 1:
        raise ValueError(f"Expected exactly one pore flag and namelist, got {flags}, {groups}")
    return text


def run_pair(solver, evidence, timeout):
    evidence.mkdir(parents=True, exist_ok=True)
    # Local fixture edits must not change the committed reference experiment.
    reference = "HEAD:TESTS/TEST_DAMBREAK/TEST_OFF_FINH/IMEX_SfloW2D.inp"
    original = subprocess.run(
        ["git", "-C", str(ROOT), "show", reference],
        check=True, capture_output=True, text=True,
    ).stdout
    variants = (("original", original), ("no_pore_pressure", without_pore_pressure(original)))
    results = []
    with tempfile.TemporaryDirectory(prefix="imex-finh-pore-") as directory:
        for name, input_text in variants:
            work = Path(directory) / name
            shutil.copytree(FIXTURE, work)
            (work / "IMEX_SfloW2D.inp").write_text(input_text, encoding="utf-8")
            subprocess.run(["sh", "./cleanFolder.sh"], cwd=work, check=True)
            output = evidence / name
            output.mkdir(parents=True, exist_ok=True)
            (output / "IMEX_SfloW2D.inp").write_text(input_text, encoding="utf-8")
            env = os.environ.copy()
            env.update(IMEX_TIMESTEP_DIAGNOSTICS="1", OMP_NUM_THREADS="1",
                       OMP_DYNAMIC="FALSE", OPENBLAS_NUM_THREADS="1")
            started = time.monotonic()
            with (output / "solver.log").open("w", encoding="utf-8") as log:
                process = subprocess.run(["timeout", "--signal=TERM", "--kill-after=15s",
                                          f"{timeout}s", str(solver)],
                                         cwd=work, env=env, stdout=log,
                                         stderr=subprocess.STDOUT, check=False)
            elapsed = time.monotonic() - started
            log_text = (output / "solver.log").read_text(encoding="utf-8", errors="replace")
            rows = extract_timesteps(log_text)
            fields = ("step", "t", "cfl_dt", "accepted_dt", "max_a_x", "max_a_y")
            with (output / "timesteps.csv").open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            result = dict(variant=name, exit_code=process.returncode,
                          elapsed_seconds=round(elapsed, 3), samples=len(rows),
                          last_step=rows[-1]["step"] if rows else None,
                          last_t=rows[-1]["t"] if rows else None,
                          min_cfl_dt=min((r["cfl_dt"] for r in rows), default=None),
                          min_accepted_dt=min((r["accepted_dt"] for r in rows), default=None),
                          max_a_x=max((r["max_a_x"] for r in rows), default=None),
                          max_a_y=max((r["max_a_y"] for r in rows), default=None),
                          invalid_timestep="INVALID_TIMESTEP" in log_text)
            results.append(result)
            print(f"{name}: status={process.returncode} samples={len(rows)} "
                  f"last_t={result['last_t']} min_cfl_dt={result['min_cfl_dt']} "
                  f"max_a_x={result['max_a_x']}", flush=True)
    with (evidence / "finh_pore_comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=int, default=120)
    args = parser.parse_args()
    solver = args.solver.resolve()
    if not solver.is_file() or not os.access(solver, os.X_OK):
        parser.error("Solver executable unavailable or not executable")
    if args.timeout_seconds <= 0:
        parser.error("--timeout-seconds must be positive")
    run_pair(solver, args.evidence.resolve(), args.timeout_seconds)


if __name__ == "__main__":
    main()
