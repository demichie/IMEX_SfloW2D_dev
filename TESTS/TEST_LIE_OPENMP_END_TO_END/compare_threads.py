#!/usr/bin/env python3
"""D20 production solver OpenMP comparison: identical fixture, 1/2/4 threads.

Runs the actual Lie timestep via two existing production examples. Rejects
missing output, nonfinite data, inconsistent output structure and unrequested
OpenMP team sizes. Never changes canonical example files.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from netCDF4 import Dataset

CASES = ("EXAMPLE_2D", "EXAMPLE_2D_PP", "EXAMPLE_2D_NEWTON", "EXAMPLE_2D_TOPO")
THREADS = (1, 2, 4)
RTOL = 5e-12
ATOL = 1e-10


def replace_one(text, expression, replacement):
    updated, count = re.subn(expression, replacement, text, count=0, flags=re.I | re.M)
    if count != 1:
        raise AssertionError(f"Expected exactly one input match for {expression}, got {count}")
    return updated


def prepare(repo, root, name):
    source_name = "EXAMPLE_2D" if name in ("EXAMPLE_2D_NEWTON", "EXAMPLE_2D_TOPO") else name
    source = repo / "EXAMPLES" / source_name
    fixture = root / "fixtures" / name
    shutil.copytree(source, fixture, ignore=shutil.ignore_patterns("OUTPUT", "__pycache__"))
    subprocess.run([sys.executable, "create_example.py", "32", "0.5", "300", "false"],
                   cwd=fixture, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                   text=True)
    input_path = fixture / "IMEX_SfloW2D.inp"
    text = input_path.read_text()
    text = replace_one(text, r"^\s*T_END\s*=.*$", " T_END=2.0E-2,")
    text = replace_one(text, r"^\s*DT_OUTPUT\s*=.*$", " DT_OUTPUT=1.0E-2,")
    # N_RK belonged to the retired Runge-Kutta integrator. The production
    # solver now uses Lie splitting and does not accept N_RK in this namelist.
    # Remove only the obsolete option from temporary test inputs.
    text = re.sub(r"^[ \t]*N_RK[ \t]*=.*$", "", text, flags=re.I | re.M)
    # The solver otherwise selects its serial path even with OMP_NUM_THREADS>1.
    if re.search(r"^\s*SERIAL_FLAG\s*=", text, re.I | re.M):
        text = replace_one(text, r"^\s*SERIAL_FLAG\s*=.*$", " SERIAL_FLAG=F,")
    else:
        text = replace_one(text, r"^\s*&RUN_PARAMETERS\b",
                           "&RUN_PARAMETERS\n SERIAL_FLAG=F,")
    text = replace_one(text, r"^\s*OUTPUT_NETCDF_FLAG\s*=.*$",
                       " OUTPUT_NETCDF_FLAG=T,")
    # Both legacy examples place exchange flags in transport namelists,
    # whereas the current solver reads them only from MASS_EXCHANGE_PARAMETERS.
    # Normalize this copied fixture, leaving the numerical settings unchanged.
    if name in CASES and not re.search(r"^\s*&MASS_EXCHANGE_PARAMETERS\b", text, re.I | re.M):
        transferred = {}
        for key in ("SETTLING_FLAG", "EROSION_COEFF", "ENTRAINMENT_FLAG"):
            pattern = rf"^\s*{key}\s*=\s*([^,\n]+)\s*,?\s*$"
            found = re.findall(pattern, text, flags=re.I | re.M)
            if len(found) != 1:
                raise AssertionError(f"{name}: expected one legacy {key}, found {len(found)}")
            transferred[key] = found[0].strip()
            text = replace_one(text, pattern, "")
        exchange = ("&MASS_EXCHANGE_PARAMETERS\n"
                    + "".join(f" {key}={value},\n" for key, value in transferred.items())
                    + " /\n")
        text = replace_one(text, r"^\s*&GAS_TRANSPORT_PARAMETERS\b",
                           exchange + "&GAS_TRANSPORT_PARAMETERS")
    if name == "EXAMPLE_2D_NEWTON":
        # Force the actual production backward-Euler Newton group.
        text = replace_one(text, r"^\s*RHEOLOGY_MODEL\s*=.*$",
                           " RHEOLOGY_MODEL='quadratic',\n RHEOLOGY_METHOD='NEWTON',\n RHEOLOGY_GROUP=2,")
    if name == "EXAMPLE_2D_TOPO":
        text = replace_one(text, r"^\s*TOPO_CHANGE_FLAG\s*=.*$",
                           " TOPO_CHANGE_FLAG=T,")
        # Dedicated active-settling fixture. A coarse solid diameter makes
        # deposition measurable against the finite-precision bed elevation;
        # disabling entrainment removes an unrelated competing mass-transfer
        # process. The ordinary EXAMPLE_2D regression remains unchanged.
        text = replace_one(text, r"^\s*DIAM_S\s*=.*$",
                           " DIAM_S=2.0E-3,")
        text = replace_one(text, r"^\s*ENTRAINMENT_FLAG\s*=.*$",
                           " ENTRAINMENT_FLAG=F,")
    input_path.write_text(text)
    return fixture


def read_fields(directory):
    paths = sorted(directory.rglob("*.nc"))
    if not paths:
        raise AssertionError(f"No NetCDF outputs under {directory}")
    fields = {}
    for path in paths:
        relative = str(path.relative_to(directory))
        with Dataset(path) as dataset:
            for name, variable in dataset.variables.items():
                if variable.dtype.kind not in "biuf":
                    continue
                raw = variable[:]
                if np.ma.isMaskedArray(raw):
                    if np.any(np.ma.getmaskarray(raw)):
                        raise AssertionError(f"Masked field {relative}/{name}")
                    raw = raw.data
                array = np.asarray(raw, dtype=np.float64)
                if not np.all(np.isfinite(array)):
                    raise AssertionError(f"Nonfinite output {relative}/{name}")
                fields[relative + "/" + name] = array.copy()
    if not fields:
        raise AssertionError("No numeric NetCDF variables")
    return fields


def run_case(solver, fixture, root, case, threads, evidence):
    target = root / "runs" / case / str(threads)
    shutil.copytree(fixture, target)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1")
    result = subprocess.run([str(solver)], cwd=target, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, timeout=240)
    log_path = evidence / f"{case}-threads-{threads}.log"
    log_path.write_text(result.stdout)
    tail = "\n".join(result.stdout.splitlines()[-35:])
    if result.returncode:
        raise AssertionError(f"{case}/{threads}: solver exited {result.returncode}; "
                             f"log={log_path}; tail:\n{tail}")
    match = re.search(r"Parallel run:\s*number of threads used\s+(\d+)", result.stdout)
    if not match:
        mode = ("non-OpenMP mode" if "Non-OpenMP simulation" in result.stdout
                else "no OpenMP team report")
        raise AssertionError(f"{case}/{threads}: {mode}; log={log_path}; "
                             f"stdout tail:\\n{tail}")
    actual = int(match.group(1))
    if actual != threads:
        raise AssertionError(f"{case}/{threads}: actual OpenMP team {actual} "
                             f"!= requested {threads}; log={log_path}; tail:\\n{tail}")
    if case == "EXAMPLE_2D_NEWTON":
        witnesses = re.findall(
            r"D20_NEWTON_WITNESS group=([0-9]+) solved_cells=([0-9]+) total_iterations=([0-9]+)",
            result.stdout)
        if not witnesses or sum(int(row[2]) for row in witnesses) <= 0:
            raise AssertionError(f"{case}/{threads}: no actual Newton iterations; log={log_path}")
    if case == "EXAMPLE_2D_TOPO":
        witnesses = re.findall(
            r"D20_BED_WITNESS group=[0-9]+ max_abs_delta=([0-9.Ee+-]+)",
            result.stdout)
        if not witnesses or max(float(value) for value in witnesses) <= 0.0:
            raise AssertionError(f"{case}/{threads}: bed did not evolve; log={log_path}")
    return read_fields(target)


def verify_full_step_rollback(solver, fixture, root, threads, evidence):
    """Require an actual production Lie rollback after a committed bed-changing group."""
    target = root / "rollback" / str(threads)
    shutil.copytree(fixture, target)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1",
               IMEX_LIE_TEST_REJECT_AFTER_GROUP="4")
    result = subprocess.run([str(solver)], cwd=target, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, timeout=240)
    log_path = evidence / f"rollback-threads-{threads}.log"
    log_path.write_text(result.stdout)
    if result.returncode == 0:
        raise AssertionError(f"rollback/{threads}: injected failure unexpectedly accepted")
    expected = ("D20_FULL_STEP_ROLLBACK_PASS", "D20_INJECTED_LATE_LIE_FAILURE")
    if not all(token in result.stdout for token in expected):
        raise AssertionError(f"rollback/{threads}: no exact full-step rollback witness; "
                             f"log={log_path}; tail:\n" +
                             "\n".join(result.stdout.splitlines()[-35:]))
    if not re.search(r"Parallel run:\s*number of threads used\s+"
                     + str(threads) + r"\b", result.stdout):
        raise AssertionError(f"rollback/{threads}: actual OpenMP team not demonstrated")
    print(f"PASS full-step injected rollback with {threads} threads")


def compare(case, baseline, candidate, threads):
    if baseline.keys() != candidate.keys():
        missing = sorted(baseline.keys() - candidate.keys())
        extra = sorted(candidate.keys() - baseline.keys())
        raise AssertionError(f"{case}/{threads}: NetCDF field mismatch: missing={missing}, extra={extra}")
    records = []
    for name, ref in baseline.items():
        other = candidate[name]
        if ref.shape != other.shape:
            raise AssertionError(f"{case}/{threads}/{name}: shape mismatch")
        difference = np.abs(other - ref)
        maximum = float(difference.max(initial=0))
        allowed = ATOL + RTOL * np.abs(ref)
        if not np.all(difference <= allowed):
            raise AssertionError(f"{case}/{threads}/{name}: max_abs={maximum:.16g}, "
                                 f"max_tolerance={float(allowed.max(initial=0)):.16g}")
        records.append({"field": name, "shape": ref.shape,
                        "max_absolute_difference": maximum})
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    args = parser.parse_args()
    repo, solver = args.repo.resolve(), args.solver.resolve()
    evidence = args.evidence.resolve()
    evidence.mkdir(parents=True, exist_ok=True)
    report = {"rtol": RTOL, "atol": ATOL, "cases": {}, "threads": THREADS}
    try:
        with tempfile.TemporaryDirectory(prefix="d20-lie-e2e-") as temp:
            root = Path(temp)
            for case in CASES:
                fixture = prepare(repo, root, case)
                baseline = None
                results = {}
                for threads in THREADS:
                    fields = run_case(solver, fixture, root, case, threads, evidence)
                    if baseline is None:
                        baseline = fields
                    else:
                        results[str(threads)] = compare(case, baseline, fields, threads)
                report["cases"][case] = {"field_count": len(baseline),
                                         "comparisons": results}
                print(f"PASS {case}: {len(baseline)} fields, teams 1/2/4")
            for threads in THREADS:
                verify_full_step_rollback(solver, root / "fixtures" / "EXAMPLE_2D_TOPO",
                                          root, threads, evidence)
            report["rollback"] = {"verified_threads": THREADS,
                                  "injected_after_group": 4}
    except Exception as error:
        report["failure"] = str(error)
        print(f"FAIL: {error}", file=sys.stderr)
        returncode = 1
    else:
        returncode = 0
    (evidence / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return returncode


if __name__ == "__main__":
    raise SystemExit(main())
