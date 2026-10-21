"""D02 cell-local characterization with independent mixture/source oracles.

A failing physical or derivative check is reported as a finding, never blessed
as an expected pass. Hazard fixtures run in separate processes with FPE traps.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np

EPS = np.finfo(float).eps
RHO_S = np.array([2500., 1800.])
CP_S = np.array([900., 1100.])
DIAM = np.array([1.e-4, 4.e-4])
SPH = np.array([.8, .6])
AIR = 101300. / (287.05 * 300.)
LAWS = (1, 2, 3, 4, 5, 6, 9, 10, 11, 12)
DEFAULTS = dict(mu=.3, xi=500., tau=.4, nu_ref=.02, visc_par=.003,
                T_ref=900., tau0=20., alpha1_coeff=1., beta1=2., alpha2=5.,
                beta2=3., Kappa=24., n_td2=.01, friction_factor=.03,
                mu_0=.2, mu_inf=.6, Fr_0=1.3, U_w=2., mu_s=.4, mu_2=.7,
                muI_inf=.005, I_0=.3, T_env=300., radiative_term_coeff=3.e-8,
                convective_term_coeff=4., hydraulic_permeability=1.e-12,
                maximum_solid_packing=.64, alpha_trans=.2, N_inh=2)


def fixture(name, model=6, h=2., temperature=500., velocity=(3., -4.),
            alpha=(.2, .1), gas=False, pressure=100., **controls):
    """Construct conservative data from prescribed volume inventories, independently."""
    a = np.array(alpha)
    carrier = 101300. / (287.05 * temperature) if gas else 1000.
    rho = a @ RHO_S[:len(a)] + (1-a.sum()) * carrier
    heat = a @ (RHO_S[:len(a)] * CP_S[:len(a)]) + (1-a.sum()) * carrier * (1005. if gas else 4180.)
    q = np.r_[rho*h, rho*h*np.array(velocity), h*heat*temperature,
              h*a*RHO_S[:len(a)], rho*h*pressure]
    config = dict(rheology_model=model, solids=len(a), gas_flag=gas,
                  bx=0., by=0., bxx=0., bxy=0., byy=0., z=0.,
                  rheology_flag=True, pore_pressure_flag=False,
                  slope_correction_flag=False, curvature_term_flag=False,
                  dynamic_permeability_flag=False, sutherland_flag=False,
                  stochastic_flag=False, f_inhibit_mode='OFF',
                  derivatives=True, operation=0, dt=.05, diagonal=1.)
    config.update(controls)
    if not config['pore_pressure_flag']:
        q=q[:-1]
    return dict(name=name, config=config, q=q)


def decode(q, c):
    """Independent mass/volume/thermal closure, including the documented dry branches."""
    n = c['solids']
    if q[0] <= EPS:
        return 0., np.zeros(2), 300., AIR, np.zeros(n), 0.
    xs = q[4:4+n]/q[0]
    cp = xs@CP_S[:n] + (1-xs.sum())*(1005. if c['gas_flag'] else 4180.)
    t = q[3]/(q[0]*cp) if q[0] > 1.e-8 else 300.
    if t <= 0:
        t = 300.
    carrier = 101300./(287.05*t) if c['gas_flag'] else 1000.
    rho = 1./(xs@(1./RHO_S[:n])+(1-xs.sum())/carrier)
    u = q[1:3]/q[0] if q[0] > 1.e-8 else math.sqrt(2)*q[0]*q[1:3]/math.sqrt(q[0]**4+1.e-32)
    return q[0]/rho, u, t, rho, xs*rho/RHO_S[:n], q[-1]/q[0]


def inhibition(a, h, c):
    mode = c['f_inhibit_mode']
    if mode == 'OFF':
        return 1.
    if mode == 'STEP':
        return float(a < c['maximum_solid_packing'])
    left = c['alpha_trans']
    if mode == 'DYNAMIC':
        left = 10**(-.28)*h**.12
        if left >= c['maximum_solid_packing']:
            left = .99*c['maximum_solid_packing']
    t = np.clip((a-left)/(c['maximum_solid_packing']-left), 0., 1.)
    # N=2 polynomial restated analytically, not copied Pascal coefficients.
    assert c['N_inh'] == 2
    return 1.-(6*t**5-15*t**4+10*t**3)


def pore_coeff(h, t, alpha, c, gas_loss=False):
    a = alpha.sum(); n = len(alpha)
    k = c['hydraulic_permeability']
    if c['dynamic_permeability_flag']:
        d = a/(alpha/DIAM[:n]).sum() if a > EPS else DIAM[:n].mean()
        weights = alpha/DIAM[:n] if gas_loss else alpha*DIAM[:n]**2
        s = weights@SPH[:n]/weights.sum() if n > 1 else SPH[0]
        k = (1-a)**3*(d*s)**2/(150*a*a)
    rho_g = 101300./(287.05*t)
    nu = 1.5e-5
    if c['gas_flag'] and c['sutherland_flag']:
        nu = 1.716e-5*(t/273.15)**1.5*(273.15+110.4)/(t+110.4)/rho_g
    beta = (1005.-287.05)/(1005.*101300.)
    return k, nu, rho_g, beta


def oracle(q, overrides):
    """Independent REAL source equations; preserves baseline conventions explicitly."""
    c = {**DEFAULTS, **overrides}
    h, u, t, rho, a, p = decode(q, c)
    g = 9.81*(rho-AIR)/rho
    vh = np.linalg.norm(u)
    v = math.sqrt(vh*vh+(u@np.array([c['bx'], c['by']]))**2) if c['slope_correction_flag'] else vh
    G = 1./(1+c['bx']**2+c['by']**2)
    C = (u[0]**2*c['bxx']+2*u[0]*u[1]*c['bxy']+u[1]**2*c['byy']) if c['curvature_term_flag'] else 0.
    implicit = np.zeros(len(q)); si = np.zeros(len(q)); model = c['rheology_model']
    if c['rheology_flag']:
        if v > 0:
            if model in (1, 12):
                xi = np.clip(c['xi']+c['z'], 1., 1.e5) if c['stochastic_flag'] else c['xi']
                implicit[1:3] = -rho*g/xi*v*u
            elif model == 2:
                implicit[1:3] = -rho*c['tau']*u/v
            elif model == 3:
                gamma = 3*c['nu_ref']/max(h, 1.e-10)*np.exp(-c['visc_par']*(t-c['T_ref']))
                implicit[1:3] = -rho*gamma*u
            elif model == 4:
                tc = t-273.15
                if tc < 20:
                    av = 1.e-3*10**(1301/(998.333+8.1855*(tc-20)+.00585*(tc-20)**2)-1.30223)
                else:
                    av = 1.002e-3*10**((1.3272*(20-tc)-.001053*(tc-20)**2)/(tc+105))
                visc = c['alpha1_coeff']*av*np.exp(c['beta1']*a.sum())
                sv = c['Kappa']*visc*v/(8*rho*9.81*max(h, 1.e-20)**2)
                std = c['n_td2']*v*v*max(h, 1.e-20)**(-4/3)
                implicit[1:3] = -9.81*rho*h*(sv+std)*u/vh
            elif model == 5:
                implicit[1:3] = -rho*1.e-3*v/(1+10*h)*u/vh
            elif model == 6:
                implicit[1:3] = -rho*c['friction_factor']*v*v*u/vh
        if model == 3:
            implicit[3] = -c['radiative_term_coeff']*(t**4-c['T_env']**4)-c['convective_term_coeff']*(t-c['T_env'])
        if h > EPS and v > 0:
            if model == 1:
                load = max(0., c['mu']*rho*(g+C)*h*math.sqrt(G))
                if c['pore_pressure_flag']:
                    load = max(0., load-c['mu']*math.sqrt(G)*p)
                si[1:3] = -max(0., load)*u/v
            elif model == 3:
                si[1:3] = -c['tau0']*u/v
            elif model == 4:
                sy = c['alpha2']*(np.exp(c['beta2']*a.sum())-1)/(9.81*rho*max(h, 1.e-20))
                si[1:3] = -9.81*rho*h*sy*u/vh
            elif model == 9 and g > 0:
                fr = v/math.sqrt(g*h)
                if c['stochastic_flag']:
                    fr = max(0., fr+(c['z'] if h > .01 else 0.))
                friction = c['mu_inf']+(c['mu_0']-c['mu_inf'])*np.exp(-fr/c['Fr_0'])
                si[1:3] = -rho*friction*h*G*g*u/vh
            elif model == 10:
                uf = v+(c['z'] if c['stochastic_flag'] else 0.)
                friction = c['mu_0'] if uf < c['U_w'] or h < .1 else c['mu_inf']+(c['mu_0']-c['mu_inf'])*c['U_w']/uf
                si[1:3] = -rho*friction*h*G*g*u/vh
            elif model in (11, 12):
                d = a.sum()/(a/DIAM[:len(a)]).sum() if a.sum() > EPS else DIAM[:len(a)].mean()
                rp = a@RHO_S[:len(a)]/a.sum() if a.sum() > EPS else RHO_S[:len(a)].mean()
                normal = max(EPS, max(0., rho*(g+C)*max(h, EPS)-(p if c['pore_pressure_flag'] else 0.))*math.sqrt(G))
                inertial = d*2.5*v/(max(h, EPS)*math.sqrt(G))/math.sqrt(normal/rp)
                friction = (c['mu_s']*c['I_0']+c['mu_2']*inertial+c['muI_inf']*inertial**2)/(c['I_0']+inertial)
                si[1:3] = -friction*normal*u/v
    if (c['pore_pressure_flag'] and h > 0 and p >= 0 and 1-a.sum() > EPS
            and (not c['dynamic_permeability_flag'] or a.sum() > EPS)):
        k, nu, rg, beta = pore_coeff(h, t, a, c)
        D = k/((1-a.sum())*nu*rg*beta)
        implicit[-1] = -rho*(math.pi/2)**2*D/max(1.e-10, h)*p*inhibition(a.sum(), h, c)
    return implicit, si


def input_text(case):
    values = []
    for key, value in case['config'].items():
        if isinstance(value, bool):
            value = '.true.' if value else '.false.'
        elif isinstance(value, str):
            value = repr(value)
        else:
            value = str(value)
        values.append(f'{key}={value}')
    values.append('input_q='+','.join(f'{x:.17e}' for x in case['q']))
    return '&fixture\n'+',\n'.join(values)+'\n/\n'


def run_probe(executable, case, directory):
    directory.mkdir(parents=True, exist_ok=True)
    data = input_text(case)
    (directory/'fixture.nml').write_text(data)
    run = subprocess.run([str(executable)], input=data, text=True, capture_output=True,
                         env={**os.environ, 'OMP_NUM_THREADS': '1', 'OMP_DYNAMIC': 'FALSE'}, timeout=20)
    (directory/'probe.log').write_text(run.stdout+run.stderr)
    parsed = {}
    for line in run.stdout.splitlines():
        tokens = line.split()
        if tokens and (tokens[0] in ('N','REAL','COMPLEX_REAL','COMPLEX_IMAG','SI','GLOBALS','STAGE','RESIDUAL','STATUS','LOSS','EXCHANGE','ODE') or tokens[0].startswith(('CS','FD','NCS'))):
            parsed[tokens[0]] = np.array([float(x) for x in tokens[1:]])
    return run, parsed


def cases():
    result = [fixture('disabled-rheology',0,rheology_flag=False)]
    for law in LAWS:
        for slope in (False, True):
            for axis, uv in [('x',(3.,0.)),('y',(0.,-4.)),('xy',(3.,-4.))]:
                result.append(fixture(f'law-{law}-slope-{int(slope)}-{axis}',law,temperature=900. if law==3 else 300.,
                    velocity=uv, slope_correction_flag=slope, bx=.35 if slope else 0., by=-.22 if slope else 0.))
        for label,h,uv in [('rest',2.,(0.,0.)),('dry',0.,(0.,0.)),('thin',1.e-12,(.01,-.02))]:
            result.append(fixture(f'law-{law}-{label}',law,h=h,velocity=uv,derivatives=False))
    for law in (1,11,12):
        for pressure in (0.,1.e3,1.e6):
            result.append(fixture(f'load-{law}-{pressure:g}',law,pressure=pressure,pore_pressure_flag=True,
                curvature_term_flag=True,bxx=.02,bxy=-.01,byy=.03,derivatives=False))
    for law in (11,12):
        result.append(fixture(f'mui-mono-{law}',law,alpha=(.3,)))
    for law in (1,12):
        for z in (-1000.,-499.,0.,99500.,100000.):
            result.append(fixture(f'xi-{law}-{z:g}',law,stochastic_flag=True,z=z))
    for law,threshold in ((9,.01),(10,.1)):
        for h in (threshold*(1-1.e-5),threshold,threshold*(1+1.e-5)):
            for z in (-100.,.2,5.):
                result.append(fixture(f'stoch-{law}-{h:g}-{z:g}',law,h=h,stochastic_flag=True,z=z,derivatives=False))
    for t in (293.1499,293.15,293.1501,280.,330.):
        result.append(fixture(f'lahar-T-{t:.12g}',4,temperature=t,derivatives=abs(t-293.15)>.01))
    for h in (.999e-10,1.e-10,1.001e-10,1.e-22,1.e-20):
        for law in (3,4):
            result.append(fixture(f'threshold-{law}-{h:g}',law,h=h,derivatives=False))
    for radiation,convection in ((0.,4.),(3.e-8,0.),(3.e-8,4.)):
        result.append(fixture(f'cooling-{radiation:g}-{convection:g}',3,temperature=900.,
            radiative_term_coeff=radiation,convective_term_coeff=convection))
    for law in (3,4,5):
        # The flat-bed implicit viscous part is differentiable at rest even
        # when the separate SI yield law is not. Test its actual derivative.
        result.append(fixture(f'viscous-rest-derivative-{law}',law,velocity=(0.,0.)))
    result.append(fixture('ambient-overwrite',3,T_env=350.))
    for speed in (1.99999,2.,2.00001):
        result.append(fixture(f'lucas-Uw-{speed:.12g}',10,h=.2,velocity=(speed,0.),derivatives=False))
    for t in (300.,900.,1100.):
        result.append(fixture(f'lava-temperature-{t:g}',3,temperature=t))
    for n in (1,2):
        for mode in ('OFF','STEP','STATIC','DYNAMIC'):
            for dynamic in (False,True):
                for suth in (False,True):
                    a = (.55,) if n==1 else (.35,.2)
                    result.append(fixture(f'pore-{n}-{mode}-{int(dynamic)}-{int(suth)}',3,gas=True,
                        temperature=700.,h=.05,alpha=a,pore_pressure_flag=True,
                        f_inhibit_mode=mode,dynamic_permeability_flag=dynamic,sutherland_flag=suth))
    for pressure in (-100.,0.,100.):
        result.append(fixture(f'pore-sign-{pressure:g}',6,gas=True,pore_pressure_flag=True,
            pressure=pressure,one_sided_pressure=pressure==0.))
    result.append(fixture('pore-only-regular',6,gas=True,rheology_flag=False,pore_pressure_flag=True))
    for a in ((0.,),(0.,0.),(1.,)):
        result.append(fixture('pore-domain-'+str(a),6,gas=True,alpha=a,pore_pressure_flag=True,
                              dynamic_permeability_flag=True,derivatives=False))
    for a in (.19999,.2,.20001,.63999,.64,.64001):
        for mode in ('STEP','STATIC'):
            result.append(fixture(f'packing-{a:g}-{mode}',6,gas=True,alpha=(a,),
                pore_pressure_flag=True,f_inhibit_mode=mode,derivatives=False))
    for law in (1,2,3,4,5,6,12):
        result.append(fixture(f'newton-{law}',law,temperature=900. if law==3 else 300.,operation=2))
    result.append(fixture('newton-pore',6,gas=True,pore_pressure_flag=True,operation=2))
    stagnation=fixture('newton-stagnation',3,h=1.e12,temperature=900.,velocity=(0.,0.),operation=2,derivatives=False)
    stagnation['expect_stagnation']=True
    result.append(stagnation)
    for law in (9,10):
        result.append(fixture(f'neutral-{law}',law,gas=True,alpha=(0.,),temperature=300.,derivatives=False))
        result.append(fixture(f'buoyant-{law}',law,gas=True,alpha=(0.,),temperature=600.,derivatives=False))
    assert len({x['name'] for x in result}) == len(result), 'fixture names must be unique'
    return result


def relative(a, b, floor=1.):
    return float(np.max(np.abs(a-b)/np.maximum(floor,np.abs(b))))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable',type=Path); parser.add_argument('output',type=Path); parser.add_argument('profile')
    args=parser.parse_args(); records=[]
    for case in cases():
        run, data = run_probe(args.executable,case,args.output/'cases'/case['name'])
        record=dict(name=case['name'],config=case['config'],q=case['q'].tolist(),exit_code=run.returncode,checks={},findings=[])
        checks=record['checks']; findings=record['findings']
        if run.returncode or 'REAL' not in data or not all(np.isfinite(x).all() for x in data.values()):
            findings.append('production source failed or did not emit a complete record')
        else:
            record['observed_globals']={k:float(v) for k,v in zip(('T_env_before','T_env_after','xi_temp'),data['GLOBALS'])}
            expected,esi=oracle(case['q'],case['config'])
            if case['name']=='packing-0.64-STEP':
                # At this discontinuity, independently rounded mass-to-volume
                # conversions can select either one-sided value. Check that
                # the output is one of those values, not an intermediate one.
                candidates=[oracle(case['q'],{**case['config'],'maximum_solid_packing':x})[0]
                            for x in (.64-1.e-12,.64+1.e-12)]
                expected=min(candidates,key=lambda x: relative(data['REAL'],x))
            checks['implicit_oracle_error']=relative(data['REAL'],expected)
            checks['si_oracle_error']=relative(data['SI'],esi)
            checks['real_complex_error']=relative(data['COMPLEX_REAL'],data['REAL'])
            checks['unperturbed_imaginary_max']=float(np.abs(data['COMPLEX_IMAG']).max())
            for key,value in checks.items():
                if not np.isfinite(value) or value > 2.e-10:
                    findings.append(key)
            if case['config']['derivatives']:
                cols=[]
                for j in range(len(case['q'])):
                    cs=data[f'CS{j+1}']; fd=[data[f'FD{j+1}_{l}'] for l in range(1,5)]
                    floor=max(1.e-12,float(np.abs(expected).max())/max(1.,abs(case['q'][j]))*1.e-5)
                    # A column norm avoids meaningless relative errors in
                    # components whose exact derivative is zero. The small
                    # absolute floor covers finite-difference cancellation.
                    errors=[float(np.max(np.abs(cs-x)))/max(floor,float(np.max(np.abs(x)))) for x in fd]
                    # Two successive step sizes must agree; no single lucky FD step suffices.
                    best=min(max(errors[k:k+2]) for k in range(3))
                    cols.append(dict(column=j+1,cs=cs.tolist(),fd_last=fd[-1].tolist(),fd_errors=errors,error=best,
                                     newton_column=(j in (1,2) or (j==3 and case['config']['rheology_model']==3)
                                                    or (j==len(case['q'])-1 and case['config']['pore_pressure_flag']))))
                record['derivatives']=cols
                checks['normalized_column_error']=max(relative(data[f'NCS{j+1}'],data[f'CS{j+1}']*max(abs(case['q'][j]),1.e-3))
                                                       for j in range(len(case['q'])))
                if checks['normalized_column_error']>2.e-10:
                    findings.append('normalized complex-step column scaling')
                if any(x['error']>3.e-5 for x in cols):
                    findings.append('complex-step versus finite-difference discrepancy')
            if case['name']=='ambient-overwrite' and data['GLOBALS'][1]!=350.:
                findings.append('configured ambient temperature overwritten')
            if case['config']['operation']==2:
                residual=float(np.max(np.abs(data['RESIDUAL'])))
                # Production tol_abs=tol_rel=1e-5, measured against the
                # initial stage residual, rather than an invented tight solve.
                bound=1.e-5+1.e-5*np.abs(case['config']['dt']*case['config']['diagonal']*expected)
                checks['stage_residual_max']=residual
                checks['stage_residual_tolerance_ratio']=float(np.max(np.abs(data['RESIDUAL'])/bound))
                if case.get('expect_stagnation'):
                    record['stagnation_detected']=bool(data['STATUS'][0]==0 and checks['stage_residual_tolerance_ratio']>1.)
                    if not record['stagnation_detected']:
                        findings.append('Newton falsely reports convergence on unrepresentable update')
                elif data['STATUS'][0]!=1 or checks['stage_residual_tolerance_ratio']>1.:
                    findings.append('Newton stage residual or convergence')
                if case['config']['rheology_model']==6 and not case['config']['pore_pressure_flag']:
                    q=case['q']; h,u,t,rho,a,p=decode(q,case['config'])
                    k=case['config']['dt']*.03/(rho*h*h); r=np.linalg.norm(q[1:3])
                    speed=2*r/(1+math.sqrt(1+4*k*r))
                    checks['quadratic_stage_oracle_error']=relative(data['STAGE'][1:3],q[1:3]*speed/r)
                    if checks['quadratic_stage_oracle_error']>1.e-5:
                        findings.append('quadratic DIRK stage oracle')
        record['status']='FAIL' if findings else 'PASS'; records.append(record)
    hazards=[fixture('pore-only',6,rheology_flag=False,pore_pressure_flag=True),
             fixture('dry-dynamic-pore',6,h=0.,pore_pressure_flag=True,dynamic_permeability_flag=True),
             fixture('zero-solids-dynamic-pore',6,alpha=(0.,0.),pore_pressure_flag=True,dynamic_permeability_flag=True),
             fixture('zero-solids-mono-pore',6,alpha=(0.,),pore_pressure_flag=True,dynamic_permeability_flag=True)]
    for case in hazards:
        case['config'].update(operation=1,derivatives=False)
        run,data=run_probe(args.executable,case,args.output/'cases'/case['name'])
        values=data.get('REAL',np.array([np.nan]))
        expected,_=oracle(case['q'],case['config'])
        ok=run.returncode==0 and np.isfinite(values).all() and relative(values,expected)<=2.e-10
        records.append(dict(name=case['name'],config=case['config'],q=case['q'].tolist(),exit_code=run.returncode,
            status='PASS' if ok else 'FAIL',findings=[] if ok else ['hazard probe: domain value or initialization failure']))
    for law in (5,6,11,12):
        for label,operation,key in [('real',1,'REAL'),('complex',7,'COMPLEX_REAL'),('si',6,'SI')]:
            case=fixture(f'rest-isolated-{law}-{label}',law,velocity=(0.,0.),operation=operation,derivatives=False)
            run,data=run_probe(args.executable,case,args.output/'cases'/case['name'])
            value=data.get(key,np.array([np.nan])); ok=run.returncode==0 and np.isfinite(value).all() and np.all(value==0.)
            records.append(dict(name=case['name'],config=case['config'],q=case['q'].tolist(),exit_code=run.returncode,
                status='PASS' if ok else 'FAIL',findings=[] if ok else ['isolated rest source path failed']))
    # Isolate Darcy gas loss: no erosion, entrainment or deposition. Check
    # volume loss and every transported conservative source independently.
    for n in (1,2):
        for mode in ('OFF','STEP','STATIC','DYNAMIC'):
            for dynamic in (False,True):
                case=fixture(f'gas-loss-{n}-{mode}-{int(dynamic)}',6,gas=True,
                    alpha=(.4,) if n==1 else (.25,.15),h=.05,temperature=700.,
                    operation=3,derivatives=False,pore_pressure_flag=True,gas_loss_flag=True,
                    f_inhibit_mode=mode,dynamic_permeability_flag=dynamic,sutherland_flag=True)
                run,data=run_probe(args.executable,case,args.output/'cases'/case['name'])
                record=dict(name=case['name'],config=case['config'],q=case['q'].tolist(),exit_code=run.returncode,
                            checks={},findings=[])
                if run.returncode or 'LOSS' not in data:
                    record['findings'].append('gas-loss production evaluation failed')
                else:
                    c={**DEFAULTS,**case['config']}; h,u,t,rho,a,p=decode(case['q'],c)
                    k,nu,rg,beta=pore_coeff(h,t,a,c,gas_loss=True)
                    requested=k/(nu*rg)/max(1.e-5,h)*math.pi/2*p*inhibition(a.sum(),h,c)
                    loss=min(requested,h*max(0.,c['maximum_solid_packing']-a.sum())/c['dt'],
                             max(0.,h*(1-a.sum())-(.4/.6)*h*a.sum())/c['dt'])
                    expected=np.zeros(len(case['q'])); expected[0]=-rg*loss
                    expected[1:3]=expected[0]*u; expected[3]=expected[0]*1005.*t
                    expected[-1]=expected[0]*p
                    record['checks']['volume_loss_error']=relative(data['LOSS'],np.array([loss]),1.e-12)
                    record['checks']['conservative_loss_error']=relative(data['EXCHANGE'],expected)
                    record['measured_loss']=float(data['LOSS'][0])
                    for key,value in record['checks'].items():
                        if not np.isfinite(value) or value>2.e-10:
                            record['findings'].append(key)
                record['status']='FAIL' if record['findings'] else 'PASS'; records.append(record)
    for label,h,uv in [('moving',2.,(3.,-4.)),('rest',2.,(0.,0.)),('dry',0.,(0.,0.))]:
        case=fixture('legacy-ode-'+label,8,h=h,velocity=uv,operation=4,derivatives=False)
        run,data=run_probe(args.executable,case,args.output/'cases'/case['name'])
        expected=case['q'].copy()
        if label=='moving':
            expected[1:3]/=1+.03*case['config']['dt']*np.linalg.norm(uv)/h
        error=relative(data.get('ODE',np.full(len(expected),np.inf)),expected)
        records.append(dict(name=case['name'],config=case['config'],q=case['q'].tolist(),exit_code=run.returncode,
            checks=dict(ode_formula_error=error),status='PASS' if run.returncode==0 and error<=2.e-10 else 'FAIL',
            findings=[] if run.returncode==0 and error<=2.e-10 else ['legacy ODE formula mismatch']))
    root=Path(__file__).resolve().parents[2]
    sources=sorted(str(p.relative_to(root)) for p in (root/'src').rglob('*') if p.suffix in ('.f90','.inc','.c'))
    sources+=sorted(str(p.relative_to(root)) for p in Path(__file__).parent.iterdir() if p.suffix in ('.py','.f90','.sh'))
    revision=subprocess.run(['git','rev-parse','HEAD'],cwd=root,capture_output=True,text=True)
    failed=[x for x in records if x['status']!='PASS']
    diagnostics=[]
    untrapped=args.executable.with_name('probe_untrapped')
    if args.profile=='optimized' and untrapped.is_file():
        for law in (5,6,11,12):
            case=fixture(f'rest-no-traps-{law}',law,velocity=(0.,0.),operation=6,derivatives=False)
            run,data=run_probe(untrapped,case,args.output/'cases'/case['name'])
            values=data.get('SI',np.array([np.nan]))
            diagnostics.append(dict(model=law,exit_code=run.returncode,finite_zero=bool(np.isfinite(values).all() and np.all(values==0.))))
    result=dict(profile=args.profile,baseline_commit=revision.stdout.strip() if revision.returncode==0 else None,
        toolchain=subprocess.check_output([os.environ.get('FC','gfortran'),'--version'],text=True).splitlines()[0],
        source_hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sources},
        untrapped_rest_diagnostics=diagnostics,
        acceptance='NOT_ACCEPTED' if failed else 'PASS',total=len(records),passed=len(records)-len(failed),failed=len(failed),
        tolerances=dict(values=2.e-10,derivatives=3.e-5,stage_residual='production 1e-5 + 1e-5*initial residual',
                        quadratic_stage_relative=1.e-5),records=records)
    (args.output/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for record in failed:
        print('FAIL:',record['name'],', '.join(record['findings']))
    print(f"D02 {args.profile}: {result['passed']}/{result['total']} fixtures pass; acceptance={result['acceptance']}",flush=True)
    return int(bool(failed))


if __name__=='__main__':
    raise SystemExit(main())
