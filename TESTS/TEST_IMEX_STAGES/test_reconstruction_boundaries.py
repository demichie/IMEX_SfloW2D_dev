"""Run production boundary reconstruction checks on 2-D and collapsed grids."""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np


def write_fixture(path, nx, ny):
    dx, dy = 0.75, 1.25
    vertices = np.zeros((ny + 1, nx + 1))
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    xx, yy = np.meshgrid(x, y)
    temperature = 300.0 + 2.0 * xx + 3.0 * yy
    state = np.zeros((ny, nx, 5))
    state[..., 0] = 1000.0
    state[..., 1] = 300.0
    state[..., 2] = -200.0
    state[..., 3] = 1000.0 * 4180.0 * temperature
    state[..., 4] = 200.0
    with path.open("w") as stream:
        stream.write(f"{nx} {ny} 0 0 5\n{dx:.17g} {dy:.17g}\n0.001\n")
        np.savetxt(stream, vertices.reshape(-1, 1), fmt="%.17g")
        np.savetxt(stream, state.reshape(-1, 1), fmt="%.17g")


def run_case(executables, root, profile, nx, ny):
    case = root / f"{profile}-{nx}x{ny}"
    case.mkdir()
    write_fixture(case / "fixture.inp", nx, ny)
    outputs = {}
    for threads in (1, 4):
        for label, executable in executables.items():
            work = case / label / f"threads-{threads}"
            work.mkdir(parents=True)
            (work / "fixture.inp").write_bytes((case / "fixture.inp").read_bytes())
            environment = {
                **os.environ,
                "OMP_NUM_THREADS": str(threads),
                "OMP_DYNAMIC": "FALSE",
            }
            run = subprocess.run(
                [str(executable), str(threads), "2", "1", "0", "0", "0", "1"],
                cwd=work,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=90,
            )
            if run.returncode:
                raise AssertionError(
                    f"{case.name}, {label}, threads={threads}:\n{run.stdout}"
                )
            if "PASS: W/E/S/N boundary modes" not in run.stdout:
                raise AssertionError(f"boundary test completion marker missing:\n{run.stdout}")
            outputs[label, threads] = (work / "boundary_traces.bin").read_bytes()
    for label in executables:
        if outputs[label, 1] != outputs[label, 4]:
            raise AssertionError(
                f"{case.name}, {label}: boundary traces differ between 1 and 4 threads"
            )
    if "baseline" in executables:
        for threads in (1, 4):
            if outputs["candidate", threads] != outputs["baseline", threads]:
                raise AssertionError(
                    f"{case.name}: candidate differs from baseline at {threads} threads"
                )
        print(
            f"PASS: {profile} baseline-identical W/E/S/N and source traces, "
            f"nx={nx}, ny={ny}, 1/4 threads",
            flush=True,
        )
    else:
        print(
            f"PASS: {profile} W/E/S/N and source reconstruction, "
            f"nx={nx}, ny={ny}, 1/4 threads",
            flush=True,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("profile")
    parser.add_argument("--baseline-executable", type=Path)
    args = parser.parse_args()
    executables = {"candidate": args.executable.resolve()}
    if args.baseline_executable:
        executables["baseline"] = args.baseline_executable.resolve()
    profile = args.profile
    with tempfile.TemporaryDirectory(prefix="imex-boundary-reconstruction-") as temp:
        root = Path(temp)
        for nx, ny in ((5, 5), (1, 5), (5, 1), (1, 1)):
            run_case(executables, root, profile, nx, ny)


if __name__ == "__main__":
    main()
