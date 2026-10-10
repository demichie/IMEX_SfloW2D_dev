"""Read-only solver diagnosis with frozen N7-B initial states (not acceptance)."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_cases as c

# Reproduce the archived diagnostic construction, not the revised acceptance
# contact. Its failed original comparison must remain observable.
c.CONTACT = c.CONTRACT['contact']

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('executable', type=Path)
parser.add_argument('profile', choices=c.CONTRACT['profiles'])
args = parser.parse_args()
executable = args.executable.resolve()
records = []
eq, p = c.CONTRACT['equilibrium'], c.P
nx, ny = eq['nx'], 1
x, y = np.meshgrid((np.arange(nx)+.5)*eq['spacing'], (np.arange(ny)+.5)*eq['spacing'])
profile = .5+.2*np.sin(.6*x)+.2*np.cos(.7*y)

for name in c.CONTRACT['closures']:
    for family in eq['families']:
        ys, base, delta = c.composition(name, profile, neutral=family=='constant-density-Q1')
        rho, _, _ = c.properties(name, ys, p['temperature'])
        if family=='constant-density-Q1':
            xv, yv = np.meshgrid(np.arange(nx+1)*eq['spacing'], np.arange(ny+1)*eq['spacing'])
            bed = 1+.1*np.cos(.4*xv)
            h = eq['eta']-.25*(bed[:-1,:-1]+bed[1:,:-1]+bed[:-1,1:]+bed[1:,1:])
            limiter = eq['density_limiter']
        else:
            bed = np.zeros((ny+1,nx+1))
            rho0, _, _ = c.properties(name, base, p['temperature'])
            h = np.sqrt((rho0-c.RHO_A)/(rho-c.RHO_A))
            limiter = eq['pressure_flat_limiter']
        initial = c.conservative(name, ys, h)
        for steps in (1,2,20):
            for n in ((2,3,4) if steps==20 else (2,)):
                label = f'{family}-{name}-RK{n}-steps{steps}'
                fields, hashes = c.launch(executable, label, name, initial, bed, n, steps, eq['dt'],
                                          eq['spacing'], eq['spacing'], (False,False), limiter)
                violations = []
                try:
                    diag = c.validate(name, fields, initial, steps)
                except AssertionError as exc:
                    violations.append(str(exc)); diag = None
                errors = np.max(np.abs(fields['q']-initial)/c.field_scales(initial),axis=(0,1))
                fractions, temp, _, _, _ = c.decode(name, fields['q'])
                max_speed = float(np.max(np.linalg.norm(fields['q'][...,1:3]/fields['q'][...,:1],axis=-1)))
                active = (fields['L'][...,1] != 0) | (fields['R'][...,1] != 0)
                # Independently calculate CU scalar numerical diffusion, separating
                # the endpoint inertial flux. Boundary/limiter modifications are
                # deliberately NOT represented as a second acceptance oracle.
                ql, qr = fields['L'], fields['R']
                fl, tl, rl, hl, _ = c.decode(name, ql)
                fr, tr, rr, hr, _ = c.decode(name, qr)
                ul, ur = ql[...,1]/ql[...,0], qr[...,1]/qr[...,0]
                cl = np.sqrt(p['gravity']*(1-c.RHO_A/rl)*hl)
                cr = np.sqrt(p['gravity']*(1-c.RHO_A/rr)*hr)
                ap = np.maximum(0,np.maximum(ul+cl,ur+cr))
                am = np.minimum(0,np.minimum(ul-cl,ur-cr))
                diffusion = ap[...,None]*am[...,None]*(qr-ql)/(ap-am)[...,None]
                record = {'case':label, 'acceptance_case':steps==20, 'threads':hashes,
                          'component_scaled_final_errors':errors.tolist(),
                          'equilibrium_gate_pass':bool(np.all(errors<=c.ROUND)),
                          'other_gate_failures':violations,'diagnostics':diag,
                          'max_velocity':max_speed,'max_temperature_change':float(np.max(np.abs(temp-300))),
                          'max_fraction_changes':np.max(np.abs(fractions-ys),axis=(0,1)).tolist(),
                          'nonzero_normal_velocity_faces':int(np.count_nonzero(active)),
                          'max_raw_CU_diffusion_nonrest':np.max(np.abs(diffusion[active]),axis=0).tolist() if np.any(active) else [0.0]*initial.shape[-1]}
                if name=='liquid' and family=='constant-density-Q1' and n==2 and steps==1:
                    advective = (ap[...,None]*ql*ul[...,None]-am[...,None]*qr*ur[...,None])/(ap-am)[...,None]
                    flux = advective+diffusion
                    component = flux[...,4:].sum(axis=-1)
                    mass = flux[...,0]
                    limited = ((mass>0)&(component>mass))|((mass<0)&(component<mass))
                    scale = np.ones_like(mass)
                    np.divide(mass,component,out=scale,where=limited)
                    flux[...,4:] *= scale[...,None]
                    rest = (ul==0)&(ur==0)
                    flux[rest,0]=0; flux[rest,3:]=0
                    # The production spatial term is +div(flux), subtracted by IMEX.
                    predicted = np.diff(flux,axis=1)/eq['spacing']
                    record['flux_diagnosis'] = {
                        'advective_component_max':np.max(np.abs(advective),axis=(0,1)).tolist(),
                        'diffusive_component_max':np.max(np.abs(diffusion),axis=(0,1)).tolist(),
                        'limiter_active_faces':np.flatnonzero(limited).tolist(),
                        'scaled_scalar_spatial_term_discrepancy':np.max(
                            np.abs(predicted[...,3:]-fields['rhs'][...,3:])/c.field_scales(initial)[3:],axis=(0,1)).tolist()
                    }
                records.append(record)
                print(label, 'errors=',errors,'u=',max_speed,'nonrest=',record['nonzero_normal_velocity_faces'],flush=True)

# Probe one prescribed contact in both signs for each closure, without claiming
# the complete acceptance inventory or its refinement gate.
s = c.CONTRACT['contact']
for name in c.CONTRACT['closures']:
    for sign in (-1,1):
        count, n = 80, 3
        dx = s['length']/count
        pulse = c.pulse_cell_means(count)
        ys, base, delta = c.composition(name,pulse)
        rho0, _, _ = c.properties(name,base,p['temperature'])
        sound = np.sqrt(p['gravity']*(1-c.RHO_A/rho0))
        coords = (np.arange(count)+.5)*dx
        rise,fall=np.clip((coords-5)/5,0,1),np.clip((35-coords)/5,0,1)
        vel=sign*(rise*rise*(3-2*rise))*(fall*fall*(3-2*fall))
        initial=c.conservative(name,ys[None],np.ones((1,count)),vel[None])
        steps=int(np.ceil(s['time']/(s['dt_dx_factor']*dx)));dt=s['time']/steps
        label=f'contact-{name}-x-sign{sign}-nx{count}-RK{n}'
        fields,hashes=c.launch(executable,label,name,initial,np.zeros((2,count+1)),n,steps,dt,dx,dx,(False,False),s['limiter'])
        failures=[]
        try: diag=c.validate(name,fields,initial,steps)
        except AssertionError as exc: failures.append(str(exc));diag=None
        expected,known=c.scalar_reference(pulse,dx,sign,sound,dt,steps,n)
        def lift(v):
            fs,_,_=c.composition(name,v[None])
            return c.conservative(name,fs,np.ones((1,count)),sign)
        expected=lift(expected); known=np.array([lift(v) for v in known])
        mask=((coords>=s['comparison_interval'][0]) & (coords<=s['comparison_interval'][1]))[None]
        errors={}
        for key in ('q','known','solved','raw_final'):
            observed=fields[key][:,mask] if fields[key].ndim==4 else fields[key][mask]
            ref=known[:,mask] if fields[key].ndim==4 else expected[mask]
            errors[key]=np.max(np.abs(observed-ref)/c.field_scales(initial),axis=tuple(range(observed.ndim-1))).tolist()
            if max(errors[key])>c.ROUND*steps:failures.append(key+': independent scalar contact mismatch')
        records.append({'case':label,'acceptance_case':True,'threads':hashes,'other_gate_failures':failures,'diagnostics':diag,'reference_errors':errors})
        print(label,'failures=',failures,'errors=',errors['q'],flush=True)
metadata = {'kind':'partial diagnostic, not N7-B closure', 'profile':args.profile,
            'revision':c.subprocess.check_output(['git','rev-parse','HEAD'],cwd=HERE,text=True).strip(),
            'test_sha256':{str(path.relative_to(HERE.parent.parent)):c.hashlib.sha256(path.read_bytes()).hexdigest()
                           for path in (HERE/'contract.json', HERE/'check_cases.py', Path(__file__),
                                        HERE.parent/'TEST_IMEX_STAGES/test_imex_stages.f90')},
            'production_sha256':{path.name:c.hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in sorted((HERE.parent.parent/'src').iterdir())
                                 if path.suffix in ('.f90','.inc') and path.is_file()},
            'records':records}
Path('diagnosis.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
print('DIAGNOSTIC ONLY: 36 paired cases recorded; this is not a passing N7-B inventory.',flush=True)


