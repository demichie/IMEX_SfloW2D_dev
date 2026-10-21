"""Source-only T2 checks through the real n_RK=2 production integrator.

Builds the existing stage observer against an up-to-date out-of-source build.
No historical Python spatial solver is imported: the central cell has zero
spatial residual in the initial uniform flat fixture for this single step.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import subprocess

import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args(); build=args.build.resolve(); out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    root=Path(__file__).resolve().parents[2]; compiler=os.environ.get('FC','gfortran')
    driver=root/'TESTS/TEST_IMEX_STAGES/test_imex_stages.f90';obj=out/'test_imex_stages.o';exe=out/'test_imex_stages'
    with (out/'build.log').open('w') as log:
        subprocess.run(['make','-C',str(build),'-j2'],stdout=log,stderr=subprocess.STDOUT,check=True)
        subprocess.run([compiler,'-O0','-g','-fcheck=all','-fopenmp','-ffpe-trap=invalid,zero,overflow',
                        '-I'+str(build/'src'),'-J'+str(out),'-c',str(driver),'-o',str(obj)],
                       cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
        text=(root/'src/Makefile.am').read_text().split('IMEX_SfloW2D_SOURCES =',1)[1].split('\n\n',1)[0]
        sources=text.replace('\\','').split();objects=[]
        for s in sources:
            if s=='IMEX_SfloW2D.f90':continue
            objects.append(str(build/'src'/str(Path(s).with_suffix('.o'))))
        libs=shlex.split(subprocess.check_output(['nf-config','--flibs'],text=True))+shlex.split(subprocess.check_output(['nc-config','--libs'],text=True))
        subprocess.run([compiler,'-fopenmp',*objects,str(obj),'-llapack',*libs,'-o',str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True)
    records=[]
    for label,law,drag,uv,t0,dt in [('quadratic',0,1,(2.,1.),300.,.1),
                                   ('coulomb-stop',1,0,(1.e-3,-2.e-3),300.,.1),
                                   ('lava-coupled',3,0,(2.,1.),900.,.1)]:
        q0=np.array([1000.,1000.*uv[0],1000.*uv[1],1000.*4180.*t0,0.])
        payloads=[]
        for threads in (1,4):
            work=out/f'{label}-t{threads}';work.mkdir(exist_ok=True)
            fields=np.tile(q0,(3,3,1))
            with (work/'fixture.inp').open('w') as f:
                f.write(f'3 3 0 0 3\n1000 1000\n{dt}\n')
                np.savetxt(f,np.zeros((16,1)),fmt='%.17g');np.savetxt(f,fields.reshape(-1,1),fmt='%.17g')
            with (work/'run.log').open('w') as log:
                run=subprocess.run([str(exe),str(threads),'2','1',str(drag),'1','0','0',str(law)],
                    cwd=work,env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_DYNAMIC':'FALSE'},stdout=log,stderr=subprocess.STDOUT)
            if run.returncode:
                records.append(dict(case=label,threads=threads,status='FAIL',reason='production integration failed',returncode=run.returncode));continue
            a=np.fromfile(work/'result.bin',dtype=np.float64);cursor=0
            def take(shape):
                nonlocal cursor
                n=int(np.prod(shape));x=a[cursor:cursor+n].reshape(shape);cursor+=n;return x
            initial=take((3,3,5)); final=take((3,3,5));take((3,3,7));take((3,3))
            known=take((2,3,3,5));solved=take((2,3,3,5));raw=take((3,3,5));statuses=take((2,3,3));stats=take((1,15))
            assert cursor==len(a) and np.isfinite(a).all()
            expected=q0.copy()
            if label=='quadratic':
                r=np.linalg.norm(q0[1:3]); k=dt/1000.
                expected[1:3]*=2./(1+math.sqrt(1+4*k*r))
            elif label=='coulomb-stop':
                expected[1:3]=0.
            else:
                lo,hi=300.,t0
                for _ in range(100):
                    temp=(lo+hi)/2
                    value=1000.*4180.*(temp-t0)+dt*(3.e-8*(temp**4-300.**4)+4*(temp-300.))
                    if value>0:hi=temp
                    else:lo=temp
                temp=(lo+hi)/2;expected[3]=1000.*4180.*temp
                speed=np.linalg.norm(q0[1:3]);si=q0[1:3]*(1-dt*20./speed)
                gamma=3*.02*math.exp(-.003*(temp-900.))
                expected[1:3]=si/(1+dt*gamma)
            def error(x,y):return float(np.max(np.abs(x-y)/np.maximum(1.,np.abs(y))))
            checks=dict(known_stage2=error(known[1,1,1],q0),source_oracle=error(raw[1,1],expected),
                        effective_increment=error(raw[1,1],solved[1,1,1]),final_projection=error(final[1,1],raw[1,1]),
                        dt_error=abs(float(stats[0,0])-dt))
            ok=(checks['source_oracle']<=1.e-5 and all(v<=1.e-10 for k,v in checks.items() if k!='source_oracle')
                and statuses[1,1,1]==1)
            records.append(dict(case=label,threads=threads,status='PASS' if ok else 'FAIL',checks=checks,
                                expected=expected.tolist(),actual=raw[1,1].tolist()))
            payloads.append(raw[1,1])
        if len(payloads)==2 and not np.array_equal(*payloads):
            records.append(dict(case=label,status='FAIL',reason='thread-dependent central state'))
    result=dict(n_RK=2,driver_sha256=hashlib.sha256(driver.read_bytes()).hexdigest(),build=str(build),records=records)
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    for r in records:print(r['status'],r['case'],r.get('threads'),r.get('checks',r.get('reason')))
    return int(any(x['status']!='PASS' for x in records))


if __name__=='__main__':
    raise SystemExit(main())
