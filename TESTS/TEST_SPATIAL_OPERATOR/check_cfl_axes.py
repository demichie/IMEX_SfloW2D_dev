"""Analytical CFL/empty-domain checks on collapsed x/y grids and a single cell."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import numpy as np
import compare_reference as spatial


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("executable",type=Path)
    cli.add_argument("profile")
    args=cli.parse_args()
    records=[]; rotated=[]
    for axis,nx,ny,dx,dy in (("x",18,1,0.75,1.25),("y",1,18,1.25,0.75),("point",1,1,0.75,1.25)):
        spatial.NX,spatial.NY,spatial.DX,spatial.DY=nx,ny,dx,dy
        q=np.zeros((ny,nx,4)); q[...,0]=700
        q[...,1]=350 if axis=="x" else 140
        q[...,2]=350 if axis=="y" else 140
        vertices=np.ones((ny+1,nx+1))
        results=[]
        for threads in (1,4):
            case=Path(f"cfl-axis-{axis}-threads-{threads}"); case.mkdir()
            spatial.write_fixture(case/"fixture.inp",vertices,q,False,False,3)
            with (case/"solver.log").open("w") as stream:
                result=subprocess.run([str(args.executable),str(threads)],cwd=case,
                    env={**os.environ,"OMP_NUM_THREADS":str(threads),"OMP_DYNAMIC":"FALSE"},
                    stdout=stream,stderr=subprocess.STDOUT,timeout=90)
            if result.returncode:
                raise AssertionError((case/"solver.log").read_text())
            steps=np.fromfile(case/"cfl_policies.bin",dtype=np.float64)
            if not np.array_equal(steps[:1],steps[1:]):
                raise AssertionError("collapsed-axis CFL policies differ")
            if axis=="point" and steps[0] != 1e6:
                raise AssertionError("no active direction failed to return MAX_DT")
            if axis!="point":
                # Liquid density is 1000; reduced gravity subtracts ambient
                # buoyancy. Uniform final traces carry exactly u_normal=0.5.
                expected=0.24*0.75/(0.5+np.sqrt(9.81*(1-1.2/1000)*0.7))
                if abs(steps[0]-expected)>8192*np.finfo(float).eps:
                    raise AssertionError("one-dimensional analytical wave-speed bound differs")
            results.append((case/"result.bin").read_bytes())
            records.append({"axis":axis,"profile":args.profile,"actual_threads":threads,
                            "dt_directional":float(steps[0]),"dt_combined":float(steps[1]),
                            "payload_sha256":hashlib.sha256(results[-1]).hexdigest()})
        if results[0]!=results[1]:
            raise AssertionError("collapsed-axis payload differs across actual 1/4 threads")
        if axis!="point": rotated.append(steps[0])
    if rotated[0] != rotated[1]:
        raise AssertionError("rotated one-dimensional timestep differs")
    Path("cfl_axes_evidence.json").write_text(json.dumps(records,indent=2)+"\n")
    print("PASS: CFL algebra on x/y/point, no propagation, caps, sentinel, actual 1/4 threads",flush=True)


if __name__=="__main__":
    main()
