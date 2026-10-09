#!/usr/bin/env python3
"""Audit N7/N8/N9 evidence at a pinned Git revision without changing numerics.

Builds and fixtures live outside the repository. A completed audit is not N9
approval: unmet requirements and failed tests remain explicit in the manifest.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from zipfile import ZipFile


TOOLS = Path(__file__).resolve().parent


def digest(path):
    value = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def probe(command):
    """Keep failed/missing version probes visible rather than inventing a version."""
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, check=False)
        return {"command": command, "returncode": result.returncode,
                "output": result.stdout.strip()}
    except OSError as error:
        return {"command": command, "error": str(error)}


def execute(command, cwd, log, environment, timeout):
    """Capture one command, terminate its process group on timeout, retain its log."""
    started = time.monotonic()
    with log.open("w") as stream:
        process = subprocess.Popen(command, cwd=cwd, env=environment,
                                   stdout=stream, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            import signal
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            code = 124
    return {"command": command, "returncode": code,
            "seconds": round(time.monotonic() - started, 3),
            "log": str(log), "log_sha256": digest(log)}


def export_revision(repo, revision, destination):
    """Export tracked sources only, without stale objects or uncommitted edits."""
    destination.mkdir()
    archive = subprocess.Popen(["git", "-C", str(repo), "archive", revision],
                               stdout=subprocess.PIPE)
    try:
        subprocess.run(["tar", "-xf", "-", "-C", str(destination)],
                       stdin=archive.stdout, check=True)
    finally:
        archive.stdout.close()
    if archive.wait() != 0:
        raise RuntimeError("git archive failed")


def run_test(source, root, executable, profile, name, threads, timeout):
    """Run an existing test; the adapter records each solver invocation before cleanup."""
    tag = f"{profile}-{name}-{threads if threads else 'embedded'}"
    case = root / "tests" / tag
    case.mkdir(parents=True)
    environment = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1",
                   "OMP_NUM_THREADS": str(threads or 1), "OMP_DYNAMIC": "FALSE"}
    command = ["sh", str(source / "TESTS" / name / "run_test.sh")]
    if executable:
        environment.update(IMEX_AUDIT_SOLVER=str(executable),
                           IMEX_AUDIT_CAPTURE=str(case / "calls"))
        if threads:
            environment["IMEX_AUDIT_FORCE_THREADS"] = str(threads)
        command.append(str(root / "adapter" / "capture_solver.sh"))
    result = execute(command, source, case / "test.log", environment, timeout)
    result.update(test=name, profile=profile, outer_threads=threads,
                  status="pass" if result["returncode"] == 0 else "fail")
    result["calls"] = []
    for path in sorted((case / "calls").glob("call-*/evidence.json")):
        call = json.loads(path.read_text())
        call["evidence_path"] = str(path)
        result["calls"].append(call)
    if executable and not result["calls"]:
        result["status"] = "fail"
        result["audit_error"] = "no solver invocation evidence"
    if any(call["problems"] for call in result["calls"]):
        result["status"] = "fail"
    print(f"{result['status'].upper()}: {tag}", flush=True)
    return result


def compare_threads(results):
    """Compare canonical q snapshots and checkpoints only for completed 1/4-thread pairs."""
    import re
    groups = {}
    for result in results:
        if result["outer_threads"] in (1, 4):
            groups.setdefault((result["profile"], result["test"]), {})[
                result["outer_threads"]] = result
    comparisons = []
    for (profile, name), pair in sorted(groups.items()):
        item = {"profile": profile, "test": name, "status": "not_completed"}
        if len(pair) == 2 and all(r["status"] == "pass" for r in pair.values()):
            def fingerprints(result):
                return [{name: record["sha256"] for name, record in call["files"].items()
                         if name == "restart.bin" or re.search(r"_\d{4}\.q_2d$", name)}
                        for call in result["calls"]]
            left, right = fingerprints(pair[1]), fingerprints(pair[4])
            has_outputs = bool(left) and all(left) and bool(right) and all(right)
            teams_match = all(
                result["calls"] and all(call["actual_threads"] == threads
                                         for call in result["calls"])
                for threads, result in pair.items())
            item["status"] = "pass" if has_outputs and teams_match and left == right else "fail"
            item["comparison"] = "bitwise SHA-256 of canonical q snapshots and restart.bin"
        comparisons.append(item)
    return comparisons


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--repo", type=Path, default=TOOLS.parents[1])
    cli.add_argument("--revision", default="HEAD")
    cli.add_argument("--output", type=Path, required=True,
                     help="new evidence directory outside the source repository")
    cli.add_argument("--netcdf", type=Path, required=True,
                     help="NetCDF Fortran prefix passed to the existing configure script")
    cli.add_argument("--jobs", type=int, default=3)
    cli.add_argument("--timeout", type=int, default=900,
                     help="per-command/test limit in seconds")
    cli.add_argument("--reference-zip", type=Path)
    cli.add_argument("--roadmap", type=Path)
    args = cli.parse_args()
    if args.jobs < 1 or args.timeout < 1:
        cli.error("jobs and timeout must be positive")
    repo, root = args.repo.resolve(), args.output.resolve()
    if root == repo or repo in root.parents:
        cli.error("evidence/build output must be outside the source repository")
    root.mkdir(parents=True, exist_ok=False)
    revision = subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", f"{args.revision}^{{commit}}"],
        text=True).strip()
    plan_path = TOOLS / "acceptance_plan.json"
    plan = json.loads(plan_path.read_text())
    source = root / "source"
    export_revision(repo, revision, source)
    adapter = root / "adapter"
    adapter.mkdir()
    for name in ("capture_solver.py", "capture_solver.sh"):
        shutil.copy2(TOOLS / name, adapter / name)
    (adapter / "capture_solver.sh").chmod(0o755)
    manifest = {
        "schema_version": 1, "started_utc": datetime.now(timezone.utc).isoformat(),
        "solver_git_revision": revision, "source_export": "git archive of tracked HEAD",
        "platform": platform.platform(), "python": sys.version,
        "audit_tools": {p.name: digest(p) for p in TOOLS.iterdir()
                        if p.is_file() and p.suffix in (".py", ".sh", ".json")},
        "source_and_fixture_sha256": {
            p.relative_to(source).as_posix(): digest(p)
            for base in (source / "src", source / "TESTS")
            for p in sorted(base.rglob("*")) if p.is_file()},
        "tools": {name: probe(command) for name, command in {
            "compiler": ["gfortran", "--version"], "autoconf": ["autoconf", "--version"],
            "automake": ["automake", "--version"],
            "netcdf_fortran": ["nf-config", "--version"],
            "netcdf_fortran_fflags": ["nf-config", "--fflags"],
            "netcdf_fortran_flibs": ["nf-config", "--flibs"],
            "netcdf_c": ["nc-config", "--all"]}.items()},
        "fixture_overrides": {"SERIAL_FLAG": False, "OMP_DYNAMIC": False},
        "unit_flag_policy": "Unit scripts compile independently with their own fixed flags; not relabeled as either production profile.",
        "profiles": {}, "tests": [], "criteria": plan["criteria"],
        "n9_decision": "open", "sign_off": "not granted by the audit runner"}
    if args.roadmap:
        manifest["roadmap"] = {"filename": args.roadmap.name,
                               "sha256": digest(args.roadmap)}
    if args.reference_zip:
        expected = {"hp_pccu_1d_core.py": "ffcb695ff8d8d2bd06e433d70b2a519c0b0012d959eb664b617deea83e8f64e2",
                    "hp_pccu_2d_core.py": "c63331366cadf9f04d802a3f058088b9c175399b2d6dc7015068a0ad110f1577"}
        with ZipFile(args.reference_zip) as archive:
            members = {Path(name).name: sha256(archive.read(name)).hexdigest()
                       for name in archive.namelist() if Path(name).name in expected}
        manifest["external_reference"] = {
            "filename": args.reference_zip.name, "sha256": digest(args.reference_zip),
            "members": members, "expected_from_gate_h_readme": expected,
            "matches_gate_h_readme": members == expected,
            "full_reference_suite_executed": False}

    def save():
        (root / "manifest.json").write_text(
            json.dumps(manifest, indent=2, allow_nan=False) + "\n")

    save()
    environment = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    bootstrap = execute(["autoreconf", "-fi"], source, root / "autoreconf.log",
                        environment, args.timeout)
    manifest["bootstrap"] = bootstrap
    executables = {}
    if bootstrap["returncode"] == 0:
        for profile, flags in plan["profiles"].items():
            print(f"BUILD: {profile}", flush=True)
            build = root / profile
            build.mkdir()
            config_env = {**environment, "FC": "gfortran", "FCFLAGS": flags}
            configured = execute([str(source / "configure"),
                                  f"--with-netcdf={args.netcdf.resolve()}"], build,
                                 build / "configure.log", config_env, args.timeout)
            record = {"flags": flags, "configure": configured}
            if configured["returncode"] == 0:
                record["make"] = execute(["make", f"-j{args.jobs}", f"FCFLAGS={flags}"],
                                         build, build / "make.log", config_env, args.timeout)
                if record["make"]["returncode"] == 0:
                    executable = build / "src" / "IMEX_SfloW2D"
                    executables[profile] = executable
                    record.update(executable_sha256=digest(executable),
                                  executable_bytes=executable.stat().st_size,
                                  linked_libraries=probe(["otool", "-L", str(executable)])
                                  if sys.platform == "darwin" else probe(["ldd", str(executable)]))
            record["status"] = "pass" if profile in executables else "fail"
            manifest["profiles"][profile] = record
            save()

    tasks = [(None, "unit_script", name, None) for name in plan["unit_tests"]]
    for profile, executable in executables.items():
        for test in plan["solver_tests"]:
            for threads in ([None] if test["threads"] == "embedded" else test["threads"]):
                tasks.append((executable, profile, test["name"], threads))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = {pool.submit(run_test, source, root, executable, profile, name,
                               threads, args.timeout): (profile, name, threads)
                   for executable, profile, name, threads in tasks}
        for future in as_completed(pending):
            try:
                result = future.result()
            except Exception as error:
                profile, name, threads = pending[future]
                result = {"test": name, "profile": profile, "outer_threads": threads,
                          "status": "fail", "audit_error": repr(error), "calls": []}
                print(f"FAIL: {profile}-{name}: {error}", flush=True)
            manifest["tests"].append(result)
            save()
    manifest["tests"].sort(key=lambda r: (r["profile"], r["test"], r["outer_threads"] or 0))
    manifest["thread_comparisons"] = compare_threads(manifest["tests"])
    by_name = {}
    for result in manifest["tests"]:
        by_name.setdefault(result["test"], []).append(result)
    for criterion in manifest["criteria"]:
        relevant = [r for name in criterion["tests"] for r in by_name.get(name, [])]
        criterion["execution_status"] = (
            "fail" if any(r["status"] == "fail" for r in relevant)
            else "pass_for_tested_scope" if relevant else "not_executed")
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    failures = [r for r in manifest["tests"] if r["status"] == "fail"]
    manifest["summary"] = {"test_runs": len(manifest["tests"]), "failed": len(failures),
                           "passed": len(manifest["tests"]) - len(failures),
                           "n9_accepted": False}
    save()
    print(json.dumps(manifest["summary"]), flush=True)
    print(f"Evidence: {root / 'manifest.json'}", flush=True)
    failed_build = bootstrap["returncode"] != 0 or len(executables) != len(plan["profiles"])
    failed_comparison = any(r["status"] == "fail" for r in manifest["thread_comparisons"])
    bad_reference = manifest.get("external_reference", {}).get("matches_gate_h_readme") is False
    return 1 if failures or failed_build or failed_comparison or bad_reference else 0


if __name__ == "__main__":
    raise SystemExit(main())
