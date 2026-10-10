"""Frozen negative controls for the exactly stationary hydrostatic correction."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import check_cases as c

GUARD = json.loads((c.HERE/'roundoff_guard_contract.json').read_text())


def run(args):
    """Require baseline equivalence for slow physical motion and real pressure forces."""
    controls = GUARD['controls']
    count, dx = controls['cells'], controls['spacing']
    x = (np.arange(count)+0.5)*dx
    profile = 0.5+0.2*np.sin(0.6*x)+0.2*np.cos(0.35)
    # Flat geometry isolates the predicate; curved/Q1 balance has its own gate.
    bed_line = np.full(count+1, controls['bed'])
    h_line = np.full(count, controls['unperturbed_depth'])
    records = []
    for name in controls['closures']:
        for axis in controls['axes']:
            transverse=controls['transverse_cells']
            ys, _, _ = c.composition(name, profile)
            ys = ys[None] if axis=='x' else ys[:,None]
            h = h_line[None] if axis=='x' else h_line[:,None]
            repeat_axis=0 if axis=='x' else 1
            ys=np.repeat(ys,transverse,axis=repeat_axis)
            h=np.repeat(h,transverse,axis=repeat_axis)
            vertices = np.broadcast_to(bed_line[None] if axis=='x' else bed_line[:,None],
                        (transverse+1,count+1) if axis=='x' else (count+1,transverse+1)).copy()
            cases = [(f'slow-{velocity}-sign{sign}', h, sign*velocity, axis)
                     for velocity in controls['velocities'] for sign in controls['directions']]
            cases.append(('tangential', h, controls['tangential_velocity'], 'y' if axis=='x' else 'x'))
            perturbation = controls['nonequilibrium_depth_perturbation']*x
            perturbed = h+(perturbation[None] if axis=='x' else perturbation[:,None])
            cases.append(('pressure-nonequilibrium', perturbed, 0.0, axis))
            for tag, depth, velocity, velocity_axis in cases:
                # BOTH axes are resolved for normal and tangential controls;
                # a one-column y fit is not supported by the pinned baseline.
                initial = c.conservative(name, ys, depth, velocity, velocity_axis)
                label = f'guard-{name}-{axis}-{tag}'
                fingerprints = {}
                for variant, executable in (('baseline',args.baseline),('candidate',args.candidate)):
                    fields, hashes = c.launch(executable, label+'-'+variant, name, initial, vertices,
                                              controls['n_RK'], controls['steps'], controls['dt'], dx, dx,
                                              (False,False), controls['limiter'])
                    fingerprints[variant] = hashes
                for payload in ('result.bin','composition.bin'):
                    baseline = Path(label+'-baseline')/'threads-1'/payload
                    candidate = Path(label+'-candidate')/'threads-1'/payload
                    if baseline.read_bytes()!=candidate.read_bytes():
                        raise AssertionError(label+': changed nonstationary/physical-force '+payload)
                speed = float(np.max(np.linalg.norm(fields['raw_final'][...,1:3]/fields['raw_final'][...,:1],axis=-1)))
                if tag=='pressure-nonequilibrium' and speed<=controls['minimum_nonequilibrium_speed']:
                    raise AssertionError(label+': physical pressure force was suppressed')
                records.append({'case':label,'fingerprints':fingerprints,'maximum_raw_final_speed':speed})
                print('PASS:',label,flush=True)
    if len(records)!=controls['expected_cases']:
        raise AssertionError('incomplete stationary-guard negative-control inventory')
    Path('roundoff_guard_evidence.json').write_text(json.dumps({
        'profile':args.profile,'contract':GUARD,
        'contract_sha256':hashlib.sha256((c.HERE/'roundoff_guard_contract.json').read_bytes()).hexdigest(),
        'cases':records,
    },indent=2,allow_nan=False)+'\n')
    print(f'PASS: {args.profile}, {len(records)} controls, baseline equality and actual 1/4 threads',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline',type=Path)
    parser.add_argument('candidate',type=Path)
    parser.add_argument('profile',choices=c.CONTRACT['profiles'])
    run(parser.parse_args())
