#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import sys

directories = sorted(
    directory for directory in os.listdir('.')
    if directory.startswith("ensemble.") and os.path.isdir(directory)
)
max_processes = 5
running = []
failures = []

if not directories:
    raise SystemExit("No ensemble.* directories found")

print(f"Start computation for {len(directories)} ensemble members")

try:
    for directory in directories:
        runner = os.path.join(directory, "run.sh")
        solver = os.path.join(directory, "IMEX_SfloW2D")
        if not (os.path.isfile(runner) and os.access(runner, os.X_OK)):
            failures.append((directory, "run.sh missing or not executable"))
            continue
        if not (os.path.isfile(solver) and os.access(solver, os.X_OK)):
            failures.append((directory, "IMEX_SfloW2D missing or not executable"))
            continue
        process = subprocess.Popen(["./run.sh"], cwd=directory)
        running.append((process, directory))
        if len(running) >= max_processes:
            process, member = running.pop(0)
            return_code = process.wait()
            if return_code:
                failures.append((member, return_code))
finally:
    for process, member in running:
        return_code = process.wait()
        if return_code:
            failures.append((member, return_code))

if failures:
    for member, return_code in failures:
        print(
            f"{member}: run.sh exited with status {return_code}",
            file=sys.stderr,
        )
    raise SystemExit(1)
