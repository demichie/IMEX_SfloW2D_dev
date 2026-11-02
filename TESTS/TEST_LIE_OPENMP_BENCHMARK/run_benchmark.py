#!/usr/bin/env python3
"""Opt-in production Lie/OpenMP throughput benchmark; independent of D20 correctness."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import random
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--cells", type=int, default=192)
    parser.add_argument("--horizon", type=float, default=0.5)
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    if args.cells < 32 or args.repetitions < 2 or args.horizon <= 0:
        parser.error("Require cells >=32, repetitions >=2, horizon >0")
    repo = args.repo.resolve()
    solver = args.solver.resolve()
    evidence = args.evidence.resolve()
    evidence.mkdir(parents=True, exist_ok=True)
    helper = repo / "TESTS/TEST_LIE_OPENMP_END_TO_END/compare_threads.py"
    spec = importlib.util.spec_from_file_location("d20_fixture", helper)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    report = {"sha": os.environ.get("GITHUB_SHA", "local"),
              "runner": os.environ.get("RUNNER_NAME", "local"),
              "cells_x": args.cells, "horizon": args.horizon,
              "repetitions": args.repetitions, "samples": {},
              "metric": "wall time including initialization and I/O"}
    try:
        with tempfile.TemporaryDirectory(prefix="lie-omp-bench-") as tmp:
            root = Path(tmp)
            fixture = module.prepare(repo, root, "EXAMPLE_2D",
                                     cells=args.cells, horizon=f"{args.horizon:.8E}")
            # Independent directories prevent overwritten outputs and warmup leakage.
            threads = (1, 2, 4, 8)
            schedule = [(n, k) for k in range(args.repetitions + 1) for n in threads]
            random.Random(20261102).shuffle(schedule)
            results = {n: [] for n in threads}
            for n, k in schedule:
                target = root / "run" / f"t{n}-r{k}"
                shutil.copytree(fixture, target)
                env = os.environ.copy()
                env.update(OMP_NUM_THREADS=str(n), OMP_DYNAMIC="FALSE",
                           OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1")
                start = time.perf_counter()
                proc = subprocess.run([str(solver)], cwd=target, env=env,
                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                      text=True, timeout=900)
                wall = time.perf_counter() - start
                log = evidence / f"threads-{n}-repeat-{k}.log"
                log.write_text(proc.stdout)
                if proc.returncode:
                    raise RuntimeError(f"Solver failed with {n} threads: {log}")
                m = re.search(r"Parallel run:\s*number of threads used\s+(\d+)", proc.stdout)
                if not m or int(m.group(1)) != n:
                    raise RuntimeError(f"Requested {n} threads not confirmed: {log}")
                if not list(target.rglob("*.nc")):
                    raise RuntimeError(f"No NetCDF output for {n} threads: {log}")
                # One untimed warmup at index 0 per thread configuration.
                if k != 0:
                    results[n].append(wall)
                print(f"threads={n}, repetition={k}, wall={wall:.4f}s", flush=True)
            baseline = statistics.median(results[1])
            rows = []
            for n in threads:
                med = statistics.median(results[n])
                speedup = baseline / med
                rows.append({"threads": n, "median_seconds": med,
                             "speedup": speedup, "efficiency": speedup / n,
                             "samples_seconds": results[n]})
            report["results"] = rows
            summary = ["# Lie OpenMP performance benchmark", "",
                       f"Commit: `{report['sha']}` · grid x={args.cells} · horizon={args.horizon}s",
                       "", "| Threads | Median wall (s) | Speedup | Efficiency |",
                       "|---:|---:|---:|---:|"]
            for row in rows:
                summary.append(f"| {row['threads']} | {row['median_seconds']:.3f} | "
                               f"{row['speedup']:.2f}x | {100*row['efficiency']:.1f}% |")
            summary += ["", "Timings include startup and I/O. GitHub-hosted runner timings "
                        "are noisy; this is an indicative comparison, not a reliable "
                        "hardware scaling certification. Verify on a dedicated host."]
            (evidence / "summary.md").write_text("\n".join(summary) + "\n")
    except Exception as exc:
        report["error"] = str(exc)
        (evidence / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        raise
    (evidence / "report.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
