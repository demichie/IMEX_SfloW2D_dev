"""Record effective inputs, actual OpenMP teams and solver output fingerprints.

This adapter never edits the repository. It sets SERIAL_FLAG=F only in the
temporary fixture created by an existing test, so requested thread counts are
really exercised. Every override and the original input are retained as evidence.
"""

from hashlib import sha256
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np


def digest(path):
    """Hash a file without loading a potentially large restart into memory."""
    value = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def fortran_float(value):
    """Read D exponents and the omitted E in formatted three-digit exponents."""
    if isinstance(value, bytes):
        value = value.decode("ascii")
    value = value.strip().replace("D", "E").replace("d", "e")
    implicit = re.fullmatch(r"([+-]?(?:\d+\.\d*|\.\d+))([+-]\d+)", value)
    if implicit:
        value = implicit.group(1) + "E" + implicit.group(2)
    return float(value)


def main():
    capture = Path(os.environ["IMEX_AUDIT_CAPTURE"])
    capture.mkdir(parents=True, exist_ok=True)
    call = capture / f"call-{len(list(capture.glob('call-*'))) + 1:03d}"
    call.mkdir()
    inputs = [p for p in Path.cwd().glob("*.inp")
              if p.name.lower() == "imex_sflow2d.inp"]
    if len(inputs) != 1:
        raise RuntimeError(f"expected one solver input, found {inputs}")
    input_path = inputs[0]
    original = input_path.read_text()
    (call / "input.original.inp").write_text(original)
    if re.search(r"(?im)^\s*SERIAL_FLAG\s*=", original):
        effective, count = re.subn(r"(?im)^\s*SERIAL_FLAG\s*=.*$",
                                  " SERIAL_FLAG=F,", original)
        if count != 1:
            raise RuntimeError("ambiguous SERIAL_FLAG declaration")
    else:
        effective, count = re.subn(r"(?i)&RUN_PARAMETERS\b",
                                  "&RUN_PARAMETERS\n SERIAL_FLAG=F,", original)
        if count != 1:
            raise RuntimeError("missing/ambiguous RUN_PARAMETERS namelist")
    overrides={"SERIAL_FLAG": False}
    selected_policy=os.environ.get("IMEX_AUDIT_CFL_POLICY")
    if selected_policy:
        # Explicit opt-in for supplementary combined-policy restart tests.
        # Never alter the default audit's fixture or hide an input override.
        if selected_policy not in ("directional_min","combined"):
            raise RuntimeError("invalid explicit audit CFL policy override")
        effective=re.sub(r"(?im)^\s*CFL_POLICY\s*=.*\n", "", effective)
        effective,count=re.subn(r"(?i)&NUMERIC_PARAMETERS\b",
            f"&NUMERIC_PARAMETERS\n CFL_POLICY='{selected_policy}',",effective)
        if count!=1:
            raise RuntimeError("missing/ambiguous NUMERIC_PARAMETERS namelist")
        overrides["CFL_POLICY"]=selected_policy
    input_path.write_text(effective)
    (call / "input.effective.inp").write_text(effective)
    environment = dict(os.environ)
    if environment.get("IMEX_AUDIT_FORCE_THREADS"):
        environment["OMP_NUM_THREADS"] = environment["IMEX_AUDIT_FORCE_THREADS"]
    requested = int(environment["OMP_NUM_THREADS"])
    command = [environment["IMEX_AUDIT_SOLVER"], *sys.argv[1:]]
    with (call / "solver.log").open("w") as log:
        process = subprocess.Popen(command, env=environment, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            log.write(line)
            sys.stdout.write(line)
        returncode = process.wait()
    log_text = (call / "solver.log").read_text()
    match = re.search(r"Parallel run: number of threads used\s+(\d+)", log_text)
    actual = int(match.group(1)) if match else None
    problems = []
    if returncode == 0 and actual != requested:
        problems.append(f"actual OpenMP team {actual} != requested {requested}")
    files = {}
    fields = {}
    for path in sorted(Path.cwd().rglob("*")):
        if path.is_file() and (path.suffix in (".inp", ".asc", ".nc", ".csv")
                               or path.name.endswith((".q_2d", ".p_2d", ".bin"))
                               or path.name == "restart.bin"):
            relative_path = str(path.relative_to(Path.cwd()))
            files[relative_path] = {"sha256": digest(path), "bytes": path.stat().st_size}
        if returncode == 0 and path.name.endswith(".q_2d"):
            try:
                values = np.loadtxt(path, ndmin=2, converters=fortran_float)
            except (ValueError, OSError) as error:
                fields[path.name] = {"parse_error": str(error)}
                problems.append(f"cannot inspect {path.name}: {error}")
                continue
            finite = bool(np.isfinite(values).all())
            fields[path.name] = {"shape": list(values.shape), "finite": finite}
            if not finite:
                problems.append(f"nonfinite conservative data in {path.name}")
            elif values.shape[1] >= 3:
                fields[path.name].update(
                    min_column_mass=float(values[:, 2].min()),
                    sum_column_mass=float(values[:, 2].sum()))
    # Effective namelist backups include initialized defaults even when the
    # user's input omits CFL_POLICY. Keep that selection in every run manifest.
    policies=[]
    for backup in Path.cwd().rglob("*.bak"):
        match=re.search(r"CFL_POLICY\s*=\s*[\"']([^\"']+)[\"']", backup.read_text(), re.I)
        if match:
            policies.append(match[1].strip())
    if len(set(policies))>1:
        problems.append("ambiguous effective CFL_POLICY backups")
    policy=policies[0] if policies else "directional_min"
    record = {"working_directory": str(Path.cwd()), "command": command,
              "returncode": returncode, "requested_threads": requested,
              "actual_threads": actual, "fixture_override": overrides,
              "input_original_sha256": sha256(original.encode()).hexdigest(),
              "input_effective_sha256": sha256(effective.encode()).hexdigest(),
              "files": files, "conservative_fields": fields, "problems": problems,
              "numerical_configuration": {"CFL_POLICY": policy,
                                           "read_from_effective_backup": bool(policies)}}
    (call / "evidence.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    if problems:
        print("AUDIT FAIL: " + "; ".join(problems), file=sys.stderr)
    return returncode if returncode > 0 else (1 if returncode or problems else 0)


if __name__ == "__main__":
    raise SystemExit(main())
