"""Production input/CFL/event contracts, without changing any numerical kernel.

Bounds fixed before execution: accepted-time algebra within 64 float64 eps;
conservative budget within 8192 eps times the initial absolute mass. Default,
sentinel and observer neutrality comparisons are byte-for-byte, not toleranced.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np

EPS = np.finfo(np.float64).eps
DT0, MAX_DT, END, OUTPUT = 0.002, 0.02, 0.095, 0.03


def check_times(values, sentinel=False):
    """Independently assemble event/growth caps from the recorded spatial proposal."""
    if values.ndim != 2 or values.shape[1] != 3 or not len(values):
        raise AssertionError("invalid accepted-time shape")
    if not np.isfinite(values).all():
        raise AssertionError("nonfinite accepted-time evidence")
    time, old, older = 0.0, DT0, DT0
    event = OUTPUT
    growth_hits = output_hits = final_hits = 0
    for t, dt, proposal in values:
        if abs(t-time) > 64*EPS or dt <= 0 or proposal <= 0 or proposal > MAX_DT+64*EPS:
            raise AssertionError("invalid contiguous positive timestep or cap")
        if sentinel and proposal != MAX_DT:
            raise AssertionError("sentinel did not return precisely MAX_DT")
        expected = min(proposal, END-t, event-t, 1.1*0.5*(old+older))
        if abs(dt-expected) > 64*EPS:
            raise AssertionError("accepted step does not obey unchanged event/growth algebra")
        growth_hits += int(dt < proposal-64*EPS and abs(dt-1.1*0.5*(old+older)) <= 64*EPS)
        output_hits += int(event < END and abs(t+dt-event) <= 64*EPS)
        final_hits += int(abs(t+dt-END) <= 64*EPS)
        time=t+dt
        older, old=old, dt
        if abs(time-event) <= 64*EPS:
            event=min(END, event+OUTPUT)
    if abs(time-END)>64*EPS or not growth_hits or output_hits != 3 or final_hits != 1:
        raise AssertionError("growth/output/final-time checks were not exercised")
    return {"accepted_steps": len(values), "growth_hits": growth_hits,
            "output_hits": output_hits, "final_hits": final_hits,
            "max_spatial_proposal": float(values[:, 2].max())}


def payloads(directory):
    """Canonical state and checkpoint bytes; exclude timestamps/diagnostic text."""
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.iterdir())
            if p.name == "restart.bin" or re.fullmatch(r"gateH_\d{4}\.q_2d", p.name)}


def prepare(directory, repo, policy, stages, cfl):
    """Generate an unchanged Q1 excavation; alter only bounded temporal settings."""
    subprocess.run(["python3", str(repo/"TESTS/TEST_PCCU_GATE_H/generate_case.py"),
                    "0.4", "g1", str(directory),
                    str(repo/"EXAMPLES/EXAMPLE_INCLINED_EXCAVATION_2D/IMEX_SfloW2D.template"),
                    "--pad-cells", "8"], check=True)
    path=directory/"IMEX_SfloW2D.inp"
    source=path.read_text()
    replacements={"T_END": END, "DT_OUTPUT": OUTPUT, "DT0": DT0,
                  "MAX_DT": MAX_DT, "N_RK": stages, "CFL": cfl,
                  "SERIAL_FLAG": "F"}
    for name, value in replacements.items():
        source, count=re.subn(rf"(?im)^\s*{name}\s*=.*$", f" {name}={value},", source)
        if count != 1:
            raise AssertionError(f"missing/ambiguous fixture setting {name}")
    source=re.sub(r"(?im)^\s*CFL_POLICY\s*=.*\n", "", source)
    if policy is not None:
        source=source.replace("&NUMERIC_PARAMETERS", f"&NUMERIC_PARAMETERS\n CFL_POLICY='{policy}',")
    path.write_text(source)


def run(executable, directory, environment, expected_error=None):
    """Require both process status and actual OpenMP evidence, or a named input error."""
    with (directory/"run.log").open("w") as log:
        result=subprocess.run([str(executable)], cwd=directory, env=environment,
                              stdout=log, stderr=subprocess.STDOUT, timeout=180)
    text=(directory/"run.log").read_text()
    if expected_error:
        if result.returncode == 0 or expected_error not in text:
            raise AssertionError(f"invalid input did not fail clearly: {directory}\n{text}")
        return
    if result.returncode:
        raise AssertionError(f"production run failed: {directory}\n{text}")
    team=re.search(r"Parallel run: number of threads used\s+(\d+)", text)
    if team is None or int(team[1]) != int(environment["OMP_NUM_THREADS"]):
        raise AssertionError("actual OpenMP team differs from requested team")


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "observer", "output", "repository"):
        cli.add_argument(name, type=Path)
    args=cli.parse_args()
    exe=args.executable.resolve(); observer=args.observer.resolve(); repo=args.repository.resolve()
    threads=int(os.environ.get("IMEX_AUDIT_FORCE_THREADS", os.environ.get("OMP_NUM_THREADS", "1")))
    env={**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE"}
    records=[]; canonical={}; proposals={}
    baseline=os.environ.get("IMEX_CFL_BASELINE")
    for stages in (2,3,4):
        for sentinel in (False, True):
            for policy in (None, "directional_min", "combined"):
                tag=f"rk{stages}-{'sentinel' if sentinel else 'adaptive'}-{policy or 'omitted'}"
                case=args.output/tag
                prepare(case, repo, policy, stages, -1 if sentinel else 0.24)
                run(exe, case, env)
                observed=args.output/(tag+"-observer")
                observed.mkdir()
                for name in ("IMEX_SfloW2D.inp", "topography_dem.asc", "gateH_0000.q_2d"):
                    shutil.copyfile(case/name, observed/name)
                run(observer, observed, env)
                fingerprint=payloads(case)
                if len(fingerprint) != 6 or fingerprint != payloads(observed):
                    raise AssertionError("observer changed canonical fields/checkpoint or output inventory")
                values=np.fromfile(observed/"accepted_times.bin", dtype=np.float64).reshape(-1,3)
                metrics=check_times(values, sentinel)
                initial=np.loadtxt(case/"gateH_0000.q_2d")[:,2]
                limit=8192*EPS*max(1.0,float(np.abs(initial).sum()))
                drift=0.0
                for path in sorted(case.glob("gateH_*.q_2d")):
                    state=np.loadtxt(path)
                    if not np.isfinite(state).all() or state[:,2].min() < -limit:
                        raise AssertionError("nonfinite or inadmissible conservative output")
                    drift=max(drift, abs(float(state[:,2].sum()-initial.sum())))
                if drift>limit:
                    raise AssertionError("closed-domain conservative mass budget failed")
                key=(stages, sentinel)
                if policy is None:
                    canonical[key]=fingerprint
                    proposals[key]=values
                elif policy == "directional_min" or sentinel:
                    if fingerprint != canonical[key] or not np.array_equal(values, proposals[key]):
                        raise AssertionError("explicit default/sentinel changed omitted-policy trajectory")
                else:
                    if not values[0,2] < proposals[key][0,2]:
                        raise AssertionError("combined policy not exercised below directional_min")
                historical=False
                if baseline and policy is None and not sentinel:
                    old=args.output/(tag+"-baseline")
                    old.mkdir()
                    for name in ("IMEX_SfloW2D.inp", "topography_dem.asc", "gateH_0000.q_2d"):
                        shutil.copyfile(case/name, old/name)
                    run(Path(baseline).resolve(),old,env)
                    if payloads(old) != fingerprint:
                        raise AssertionError("omitted policy changed the pinned historical executable")
                    historical=True
                backup=(case/"gateH.bak").read_text()
                selected=policy or "directional_min"
                if not re.search(rf"CFL_POLICY\s*=\s*[\"']{selected}\s*[\"']",backup,re.I):
                    raise AssertionError("effective input backup lacks selected CFL_POLICY")
                records.append({"case": tag, "actual_threads": threads, "policy": selected,
                                "parameter_omitted": policy is None, "n_RK": stages,
                                "sentinel": sentinel, "observer_byte_identical": True,
                                "historical_baseline_byte_identical": historical,
                                "canonical_sha256": fingerprint,
                                "accepted_times_sha256": hashlib.sha256((observed/"accepted_times.bin").read_bytes()).hexdigest(),
                                "input_sha256": hashlib.sha256((case/"IMEX_SfloW2D.inp").read_bytes()).hexdigest(),
                                "mass_drift": drift, "mass_limit": limit, **metrics})
                print(f"PASS: {tag}, actual {threads} threads, proposal/event/growth/output/checkpoint",flush=True)
    for tag,policy,cfl,message in (("unknown","invalid",0.24,"CFL_POLICY must be"),
                                   ("negative","combined",-2,"CFL value outside"),
                                   ("above_bound","combined",0.26,"CFL value outside")):
        case=args.output/tag
        prepare(case,repo,policy,2,cfl)
        # Expected input failures are negative controls, not successful solver
        # trajectories: do not send them to the audit's success-only adapter.
        run(Path(env.get("IMEX_AUDIT_SOLVER",str(exe))).resolve(),case,env,message)
    evidence={"status":"passed", "actual_threads":threads, "cases":records,
              "negative_controls":["unknown policy","CFL=-2","CFL=0.26 in 2D"],
              "baseline_executable_sha256": hashlib.sha256(Path(baseline).read_bytes()).hexdigest() if baseline else None,
              "time_bound":64*EPS, "budget_factor":8192,
              "n9_approved":False}
    (args.output/"evidence.json").write_text(json.dumps(evidence,indent=2,allow_nan=False)+"\n")
    print("PASS: production CFL policies and accepted-time contract",flush=True)


if __name__ == "__main__":
    main()
