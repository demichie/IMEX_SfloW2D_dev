"""Compare default accepted/proposed steps with the frozen pre-CFL executable.

Complement check_cases.py's uninstrumented historical field comparison. The
baseline observer must independently preserve its own canonical outputs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import numpy as np
from check_cases import payloads, run


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    for name in ("case_root","baseline","repository"):
        cli.add_argument(name,type=Path)
    args=cli.parse_args()
    baseline=args.baseline.resolve(); repo=args.repository.resolve(); root=args.case_root.resolve()
    evidence=json.loads((root/"evidence.json").read_text())
    env={**os.environ,"OMP_NUM_THREADS":str(evidence["actual_threads"]),"OMP_DYNAMIC":"FALSE"}
    env.pop("IMEX_AUDIT_SOLVER",None)
    build=root/"historical-time-observer-build"
    subprocess.run(["python3",str(repo/"TESTS/TEST_PCCU_GATE_H/build_dt_observer.py"),
                    str(baseline),str(build),str(repo),"--include-spatial"],env=env,check=True)
    records=[]
    for stages in (2,):
        tag=f"rk{stages}-adaptive-omitted"
        original=root/(tag+"-baseline")
        observed=root/(tag+"-baseline-observer")
        observed.mkdir()
        for name in ("IMEX_SfloW2D.inp","topography_dem.asc","gateH_0000.q_2d"):
            shutil.copyfile(original/name,observed/name)
        run(build/"observe_accepted_times",observed,env)
        if payloads(observed)!=payloads(original):
            raise AssertionError("historical observer changed canonical fields/checkpoint")
        left=(root/(tag+"-observer")/"accepted_times.bin").read_bytes()
        right=(observed/"accepted_times.bin").read_bytes()
        if left!=right:
            raise AssertionError("default accepted/proposed times differ from pre-CFL implementation")
        records.append({"n_RK":stages,"actual_threads":evidence["actual_threads"],
                        "historical_observer_byte_identical":True,
                        "accepted_and_proposed_steps_byte_identical":True,
                        "accepted_steps":len(np.frombuffer(left,dtype=np.float64))//3,
                        "accepted_times_sha256":hashlib.sha256(left).hexdigest()})
    result={"status":"passed","baseline_executable_sha256":hashlib.sha256(baseline.read_bytes()).hexdigest(),
            "records":records}
    (root/"baseline_times_evidence.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS: default proposed/accepted steps and canonical fields match 691a3e21 byte for byte",flush=True)


if __name__=="__main__":
    main()
