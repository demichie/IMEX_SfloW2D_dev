"""Observe accepted steps in a separate run; reject any change to canonical fields."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


def capture(executable, case):
    """Copy only initial inputs, then verify the I/O hook is numerically inert."""
    directory = case/"accepted-time-observer"
    directory.mkdir()
    for name in ("IMEX_SfloW2D.inp", "topography_dem.asc", "gateH_0000.q_2d"):
        shutil.copyfile(case/name, directory/name)
    with (directory/"run.log").open("w") as stream:
        subprocess.run([str(executable.resolve())], cwd=directory,
                       env=os.environ.copy(), stdout=stream, stderr=subprocess.STDOUT,
                       check=True, timeout=600)
    original = (case/"gateH_0001.q_2d").read_bytes()
    observed = (directory/"gateH_0001.q_2d").read_bytes()
    if original != observed:
        raise AssertionError("accepted-time observer changed canonical conservative fields")
    match = re.search(r"Parallel run: number of threads used\s+(\d+)",
                      (directory/"run.log").read_text())
    requested = int(os.environ.get("OMP_NUM_THREADS", "1"))
    actual = int(match[1]) if match else None
    if actual != requested:
        raise AssertionError(f"observational OpenMP team {actual} != requested {requested}")
    shutil.copyfile(directory/"accepted_times.bin", case/"accepted_times.bin")
    evidence = {"canonical_output_byte_identical": True, "actual_threads": actual,
                "canonical_output_sha256": hashlib.sha256(original).hexdigest(),
                "accepted_times_sha256": hashlib.sha256((case/"accepted_times.bin").read_bytes()).hexdigest(),
                "numerical_modules_instrumented": False}
    (case/"accepted_time_evidence.json").write_text(json.dumps(evidence, indent=2)+"\n")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("executable", type=Path)
    cli.add_argument("case", type=Path)
    args = cli.parse_args()
    capture(args.executable, args.case)
