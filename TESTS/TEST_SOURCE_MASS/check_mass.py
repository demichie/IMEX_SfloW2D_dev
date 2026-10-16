#!/usr/bin/env python3
"""Compare the mass in the domain with the mass the source must have injected.

For each output of a case (the .q_2d files), the solid mass, or for the
lateral radial source the mixture mass, is summed over the domain and compared
with the time integral of the prescribed source. The integrals are evaluated
here independently of the solver: each window is split at its corners and
every piece is integrated with a 16-point Gauss-Legendre rule, exact for the
linear and cosine ramps to round-off.

Usage: check_mass.py CASE_DIR [TOLERANCE]
"""
import glob
import json
import math
import os
import sys

import numpy as np

ALPHA_S, RHO_S = 0.02, 2000.0
GL_X, GL_W = np.polynomial.legendre.leggauss(16)


def window(s, per, act, ramp, phase, cosine):
    """Point value of a periodic source window at time s."""
    u = math.fmod(s + phase, per)
    if ramp <= 0.0:
        return 1.0 if u <= act else 0.0
    if u < ramp:
        return 0.5 * (1.0 - math.cos(math.pi * u / ramp)) if cosine else u / ramp
    if u <= act - ramp:
        return 1.0
    if u <= act:
        if cosine:
            return 0.5 * (1.0 + math.cos(math.pi * (u - act + ramp) / ramp))
        return (act - u) / ramp
    return 0.0


def integral(f, corners, t):
    """Integral of f over [0, t], split at the given corner times."""
    pts = sorted({0.0, t} | {c for c in corners if 0.0 < c < t})
    total = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        h = 0.5 * (b - a)
        total += h * sum(w * f(a + h * (1.0 + x)) for x, w in zip(GL_X, GL_W))
    return total


def window_integral(t, per, act, ramp, phase, cosine, cutoff=None):
    end = t if cutoff is None else min(t, cutoff)
    if end <= 0.0:
        return 0.0
    corners = []
    k = -2
    while k * per - phase <= end + per:
        for c in (0.0, ramp, act - ramp, act):
            corners.append(k * per - phase + c)
        k += 1
    return integral(lambda s: window(s, per, act, ramp, phase, cosine), corners, end)


def expected(meta, t):
    kind = meta['kind']
    if kind == 'fissures':
        total = 0.0
        for q in meta['fissures']:
            rate = ALPHA_S * RHO_S * q['vel'] * q['length'] * q['width']
            if q['tp'] is None:
                total += rate * t
            else:
                per, act, ramp, phase = q['tp']
                total += rate * window_integral(t, per, act, ramp, phase, False)
        return total
    if kind == 'lateral':
        per, act, ramp, cutoff = meta['tp']
        return meta['mfr'] * window_integral(t, per, act, ramp, 0.0, True, cutoff)
    if kind == 'bottom':
        per, act, ramp, phase = meta['tp']
        return (ALPHA_S * RHO_S * meta['vel'] * meta['area']
                * window_integral(t, per, act, ramp, phase, False))
    if kind == 'intervals':
        ts, vs = meta['t'], meta['vel']
        total = 0.0
        for i, (a, v) in enumerate(zip(ts, vs)):
            b = ts[i + 1] if i + 1 < len(ts) else float('inf')
            total += v * max(0.0, min(t, b) - max(0.0, a))
        return ALPHA_S * RHO_S * meta['area'] * total
    raise ValueError(kind)


def main():
    case = sys.argv[1]
    tol = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0e-9
    meta = json.load(open(os.path.join(case, 'sources.json')))
    files = sorted(glob.glob(os.path.join(case, meta['name'] + '_[0-9][0-9][0-9][0-9].q_2d')))
    if len(files) < 2:
        print('FAIL: %s: no output' % meta['name'])
        return 1
    # q_2d columns: x, y, then the conservative variables; with one solid
    # class and no additional gas: rho*h, rho*h*u, rho*h*v, energy, solid
    column = 2 if meta['kind'] == 'lateral' else 6
    area = meta['dx'] ** 2
    worst, failed = 0.0, False
    final = expected(meta, meta['t_end'])
    for i, f in enumerate(files[1:], 1):
        t = i * meta['dt_out']
        q = np.loadtxt(f)
        mass = q[:, column].sum() * area
        e = expected(meta, t)
        if e > 1.0e-9 * final:
            err = mass / e - 1.0
            bad = abs(err) > tol
            print('  %-24s t = %5.1f s  expected %.10e kg  in domain %.10e kg  %+.2e%s'
                  % (meta['name'], t, e, mass, err, '  FAIL' if bad else ''))
            worst = max(worst, abs(err))
        else:
            bad = abs(mass) > 1.0e-9 * final
            print('  %-24s t = %5.1f s  nothing injected yet, in domain %.3e kg%s'
                  % (meta['name'], t, mass, '  FAIL' if bad else ''))
        failed = failed or bad
    print('%s: %s (worst relative error %.2e, tolerance %.0e)'
          % ('FAIL' if failed else 'PASS', meta['name'], worst, tol))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
