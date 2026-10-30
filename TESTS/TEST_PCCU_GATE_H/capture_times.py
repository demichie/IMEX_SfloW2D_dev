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
    def snapshot(run_directory):
        snapshots = run_directory/"OUTPUT/gateH/snapshots"
        matches = sorted(snapshots.glob("*/gateH_0001.q_2d"))
        if len(matches) != 1:
            manifest = run_directory/"OUTPUT/gateH/metadata/output_manifest.csv"
            manifest_text = manifest.read_text() if manifest.exists() else "<missing>"
            logfile = run_directory/"run.log"
            log_tail = "\\n".join(logfile.read_text(errors="replace").splitlines()[-50:]) if logfile.exists() else "<missing>"
            found = sorted(str(path.relative_to(run_directory))
                           for path in (run_directory/"OUTPUT").rglob("*.q_2d"))
            raise AssertionError(
                f"expected one gateH_0001.q_2d, found {len(matches)} in {run_directory}\\n"
                f"available snapshots: {found}\\n"
                f"output manifest:\\n{manifest_text}\\n"
                f"solver log tail:\\n{log_tail}")
        return matches[0]

    original_path = snapshot(case)
    observed_path = snapshot(directory)
    original = original_path.read_bytes()
    observed = observed_path.read_bytes()
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
