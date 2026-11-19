#!/usr/bin/env python3
"""Compare early CFL/timestep behavior of the four unchanged FINH dam-break fixtures.

This diagnostic intentionally allows a short per-case wall-clock timeout:
timing out is evidence, not a reason to change the physical simulation horizon.
No GitHub Action is dispatched by this script.
"""

import argparse
import csv
import json
import os
from pathlib import Path
import re
import subprocess
import time


ROOT = Path(__file__).resolve().parents[1]
CASES = ("TEST_OFF_FINH", "TEST_STEP_FINH", "TEST_STATIC_FINH", "TEST_DYNAMIC_FINH")
WITNESS = re.compile(
    r"^TIMESTEP_DIAG step=(\d+)\s+t=\s*(\S+)\s+"
    r"cfl_dt=\s*(\S+)\s+accepted_dt=\s*(\S+)\s+"
    r"max_a_x=\s*(\S+)\s+max_a_y=\s*(\S+)\s*$",
    re.MULTILINE,
)


def extract_timesteps(log):
    """Read actual scientific-notation quantities; never use rounded 'dt = 0.00000'."""
    rows = []
    for match in WITNESS.finditer(log):
        step = int(match.group(1))
        values = [float(value.replace("D", "E")) for value in match.groups()[1:]]
        rows.append(dict(zip(
            ("step", "t", "cfl_dt", "accepted_dt", "max_a_x", "max_a_y"),
            (step, *values),
        )))
    return rows


def run_comparison(solver, evidence, timeout_seconds, selected):
    evidence.mkdir(parents=True, exist_ok=True)
    report = []
    for case in selected:
        directory = evidence / case
        directory.mkdir(parents=True, exist_ok=True)
        environment = os.environ.copy()
        environment.update(
            IMEX_TIMESTEP_DIAGNOSTICS="1",
            IMEX_TEST_LOG_DIR=str(directory),
            IMEX_TEST_CASE_TIMEOUT=f"{timeout_seconds}s",
            OMP_NUM_THREADS="1",
            OMP_DYNAMIC="FALSE",
            OPENBLAS_NUM_THREADS="1",
        )
        started = time.monotonic()
        try:
            result = subprocess.run(
                ["sh", str(ROOT / "TESTS" / "run_tests.sh"), str(solver), case],
                cwd=ROOT,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=timeout_seconds + 45,
                check=False,
            )
            status = result.returncode
            output = result.stdout
        except subprocess.TimeoutExpired as exc:
            status = 124
            output = (exc.stdout or b"")
            if isinstance(output, bytes):
                output = output.decode("utf-8", "replace")
            output += "\nDIAGNOSTIC_WRAPPER_TIMEOUT\n"
        elapsed = time.monotonic() - started
        (directory / "runner.log").write_text(output, encoding="utf-8")
        solver_log = directory / f"{case}.solver.log"
        content = solver_log.read_text(encoding="utf-8", errors="replace") if solver_log.exists() else ""
        timesteps = extract_timesteps(content)
        row = {
            "case": case,
            "exit_code": status,
            "elapsed_seconds": round(elapsed, 3),
            "samples": len(timesteps),
            "last_step": timesteps[-1]["step"] if timesteps else None,
            "last_t": timesteps[-1]["t"] if timesteps else None,
            "min_cfl_dt": min((step["cfl_dt"] for step in timesteps), default=None),
            "min_accepted_dt": min((step["accepted_dt"] for step in timesteps), default=None),
            "max_a_x": max((step["max_a_x"] for step in timesteps), default=None),
            "max_a_y": max((step["max_a_y"] for step in timesteps), default=None),
            "invalid_timestep": "INVALID_TIMESTEP" in content,
        }
        report.append(row)
        with (directory / "timesteps.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=("step", "t", "cfl_dt", "accepted_dt", "max_a_x", "max_a_y")
            )
            writer.writeheader()
            writer.writerows(timesteps)
        print(
            f"{case}: status={status} elapsed={elapsed:.1f}s "
            f"samples={len(timesteps)} last_t={row['last_t']} "
            f"min_cfl_dt={row['min_cfl_dt']} min_accepted_dt={row['min_accepted_dt']}",
            flush=True,
        )
    (evidence / "finh_comparison.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    with (evidence / "finh_comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(report[0]))
        writer.writeheader()
        writer.writerows(report)
    if any(row["samples"] == 0 for row in report):
        raise RuntimeError("Missing TIMESTEP_DIAG samples: inspect per-case runner and solver logs")
    print("Diagnostic comparison recorded; timeouts are reported, not treated as numerical acceptance.")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, required=True,
                        help="Existing compiled solver, optimized or strict-debug")
    parser.add_argument("--evidence", type=Path, required=True,
                        help="Directory for solver logs and numerical comparison")
    parser.add_argument("--timeout-seconds", type=int, default=90,
                        help="Maximum wall time per solver run (default: 90 seconds)")
    parser.add_argument("--case", action="append", choices=CASES,
                        help="Choose one or more modes; default: all four")
    args = parser.parse_args()
    if args.timeout_seconds <= 0:
        parser.error("--timeout-seconds must be positive")
    solver = args.solver.resolve()
    if not solver.is_file() or not os.access(solver, os.X_OK):
        parser.error(f"Solver executable unavailable: {solver}")
    run_comparison(solver, args.evidence.resolve(), args.timeout_seconds, args.case or CASES)


if __name__ == "__main__":
    main()
