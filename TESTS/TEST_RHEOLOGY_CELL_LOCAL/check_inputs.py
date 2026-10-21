"""D02 parser-domain probes through the actual solver; plain STOP can exit zero."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args();exe=args.executable.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    root=Path(__file__).resolve().parents[2]
    subprocess.run(['python3','-B',str(root/'TESTS/TEST_SOURCE_MASS/make_cases.py'),str(out/'base')],check=True)
    base=out/'base/fissures';template=(base/'IMEX_SfloW2D.inp').read_text()
    # These fixed-name directories contain only generated test artifacts;
    # make the documented command repeatable without deleting user paths.
    template=template.replace('&RUN_PARAMETERS\n','&RUN_PARAMETERS\n OUTPUT_OVERWRITE=T,\n')
    template=template.replace('RHEOLOGY_FLAG = F','RHEOLOGY_FLAG = T')
    template=re.sub(r'T_END=\s*40\.0D0','T_END= 0.0001D0',template)
    cases=[('valid-voellmy',1,'mu=.3, xi=500',False),
           ('unsupported-7',7,'',True),('disabled-active-0',0,'',True),
           ('legacy-8',8,'friction_factor=.03',False),
           ('missing-xi',1,'mu=.3',True),('missing-mu0',9,'mu_inf=.6, Fr_0=1.3',True),
           ('zero-Fr0',9,'mu_0=.2, mu_inf=.6, Fr_0=0',True),
           ('zero-Uw',10,'mu_0=.6, mu_inf=.2, U_w=0',True),
           ('missing-muinf',10,'mu_0=.6, U_w=2',True),
           ('all-negative-froude',9,'mu_0=-1, mu_inf=-1, Fr_0=-1',True),
           ('zero-coulomb',1,'mu=0, xi=500',False),
           ('zero-froude-friction',9,'mu_0=0, mu_inf=0, Fr_0=1.3',False),
           ('zero-lucas-friction',10,'mu_0=0, mu_inf=0, U_w=2',False),
           ('nan-mu',1,'mu=NaN, xi=500',True),
           ('infinite-xi',1,'mu=.3, xi=Infinity',True),
           ('nan-Fr0',9,'mu_0=.2, mu_inf=.6, Fr_0=NaN',True),
           ('infinite-Uw',10,'mu_0=.6, mu_inf=.2, U_w=Infinity',True)]
    cases.extend([
        ('nan-plastic',2,'tau=NaN',True),
        ('nan-lava',3,'nu_ref=NaN, tau0=20, visc_par=.003, T_ref=900',True),
        ('nan-lahar',4,'alpha1_ref=NaN',True),
        ('infinite-quadratic',6,'friction_factor=Infinity',True),
        ('nan-mui',11,'mu_s=.4, mu_2=.7, I_0=.3, muI_inf=NaN',True),
        ('nan-mui-voellmy',12,'mu_s=.4, mu_2=.7, I_0=.3, muI_inf=.005, xi=NaN',True)])
    for name in ('ambient-default','ambient-set','ambient-zero','ambient-nan','ambient-inf'):
        cases.append((name,3,'nu_ref=.02, tau0=20, visc_par=.003, T_ref=900',
                      name in ('ambient-zero','ambient-nan','ambient-inf')))
    records=[]
    for name,model,params,should_reject in cases:
        work=out/name;work.mkdir(exist_ok=True)
        shutil.copy2(base/'flat_dem.asc',work/'flat_dem.asc')
        text=template
        if model==3:
            text=text.replace("'froude_blend'","'h_x'")
            ambient={'ambient-default':'','ambient-set':'T_ENV=350,',
                     'ambient-zero':'T_ENV=0,','ambient-nan':'T_ENV=NaN,',
                     'ambient-inf':'T_ENV=Infinity,'}.get(name,'')
            text=text.replace('&TEMPERATURE_PARAMETERS\n /',
                '&TEMPERATURE_PARAMETERS\n EMISSIVITY=.8, ATM_HEAT_TRANSF_COEFF=4, EXP_AREA_FRACT=.5, '+ambient+'\n /')
        (work/'IMEX_SfloW2D.inp').write_text(text+f'\n&RHEOLOGY_PARAMETERS\n rheology_model={model}, {params}\n/\n')
        with (work/'run.log').open('w') as f:
            run=subprocess.run([str(exe)],cwd=work,stdout=f,stderr=subprocess.STDOUT,
                env={**os.environ,'OMP_NUM_THREADS':'1','OMP_DYNAMIC':'FALSE'},timeout=30)
        log=(work/'run.log').read_text()
        token=('ERROR: invalid T_ENV in TEMPERATURE_PARAMETERS' if name.startswith('ambient-')
               else 'ERROR: problem with namelist RHEOLOGY_PARAMETERS')
        rejected=token in log
        # A process error elsewhere is not evidence that rheology validation rejected input.
        reached_end='Total elapsed real time' in log
        passed=(rejected if should_reject else (not rejected and reached_end and run.returncode==0))
        ambient_temperature=None
        if name in ('ambient-default','ambient-set') and reached_end:
            match=re.search(r'Cooling ambient temperature \(K\)\s*=\s*([^,\s]+)',log)
            ambient_temperature=float(match[1].replace('D','E')) if match else None
            passed=passed and ambient_temperature==(350. if name=='ambient-set' else 300.)
        records.append(dict(name=name,model=model,should_reject=should_reject,rejected=rejected,
                            reached_end=reached_end,exit_code=run.returncode,ambient_temperature=ambient_temperature,
                            status='PASS' if passed else 'FAIL'))
    result=dict(records=records,acceptance='PASS' if all(r['status']=='PASS' for r in records) else 'NOT_ACCEPTED')
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    for r in records:print(r['status'],r['name'],'rejected=',r['rejected'],'exit=',r['exit_code'])
    return int(result['acceptance']!='PASS')


if __name__=='__main__':
    raise SystemExit(main())
