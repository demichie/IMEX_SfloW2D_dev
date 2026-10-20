#!/usr/bin/env python3
"""Write the TEST_SOURCE_MASS cases: one input file and one flat DEM each.

Every case runs on a flat bed with settling, erosion and entrainment off, so
the mass in the domain at an output time is the mass the source injected.
The source parameters are kept in sources.json for check_mass.py.

Usage: make_cases.py WORK_DIR
"""
import json
import math
import os
import sys

TEMPLATE = """&RUN_PARAMETERS
 RUN_NAME="{name}",
 RESTART=F,
 T_START= 0.0D0,
 T_END= {t_end:.1f}D0,
 DT_OUTPUT= {dt_out:.1f}D0,
 OUTPUT_CONS_FLAG=T,
 OUTPUT_ESRI_FLAG=F,
 OUTPUT_NETCDF_FLAG=F,
 OUTPUT_RUNOUT_FLAG=F,
 VERBOSE_LEVEL=0,
 /
&NEWRUN_PARAMETERS
 TOPOGRAPHY_FILE = 'flat_dem.asc',
 X0= {x0:.1f}D0,
 Y0= {y0:.1f}D0,
 COMP_CELLS_X= {n},
 COMP_CELLS_Y= {n},
 CELL_SIZE = {dx:.1f}D0,
 N_SOLID = 1,
 N_ADD_GAS = 0,
 RHEOLOGY_FLAG = F,
 TOPO_CHANGE_FLAG = F,
{flags} COLLAPSING_VOLUME_FLAG = F,
 /
{bcs}&TEMPERATURE_PARAMETERS
 /
&NUMERIC_PARAMETERS
 KT_TRACE_RECONSTRUCTION = 'froude_blend',
 DT0 = 1.0E-4,
 MAX_DT= 0.5D0,
 CFL= 0.2D0,
 LIMITER= 7*3,
 THETA= 1.3D0,
 RECONSTR_COEFF= 1.0D0,
 N_RK= 2,
 /
&EXPL_TERMS_PARAMETERS
 GRAV= 9.81D0,
 /
&SOLID_TRANSPORT_PARAMETERS
 RHO_S = 2000.0D0,
 DIAM_S = 1.0E-4,
 SP_HEAT_S = 1000.0D0,
 SETTLING_FLAG = F,
 EROSION_COEFF = 0.0D0,
 /
{source}&GAS_TRANSPORT_PARAMETERS
 SP_HEAT_A = 998.0D0,
 SP_GAS_CONST_A = 287.051D0,
 KIN_VISC_A = 1.48D-5,
 PRES = 101300.0D0,
 T_AMBIENT = 300.0D0,
 ENTRAINMENT_FLAG = F,
 /
"""


def bcs():
    out = ''
    for side, c in (('WEST', 'W'), ('EAST', 'E'), ('SOUTH', 'S'), ('NORTH', 'N')):
        out += '&%s_BOUNDARY_CONDITIONS\n' % side
        for v in ('H', 'HU', 'HV'):
            out += ' %s_BC%s%%FLAG=1, %s_BC%s%%VALUE=0.0D0,\n' % (v, c, v, c)
        out += ' XS_BC%s(1)%%FLAG=1, XS_BC%s(1)%%VALUE=0.0D0,\n' % (c, c)
        out += ' T_BC%s%%FLAG=1, T_BC%s%%VALUE=0.0D0,\n /\n' % (c, c)
    return out


def flat_dem(path, x0, y0, n, dx):
    # two cells of margin around the computational domain
    m = n + 4
    with open(path, 'w') as f:
        f.write('ncols %d\nnrows %d\nxllcorner %.1f\nyllcorner %.1f\ncellsize %.1f\n'
                'NODATA_value -9999\n' % (m, m, x0 - 2 * dx, y0 - 2 * dx, dx))
        row = ' '.join(['0.0'] * m) + '\n'
        for _ in range(m):
            f.write(row)


def fissure(xc, yc, length, ang_deg):
    dx = 0.5 * length * math.cos(math.radians(ang_deg))
    dy = 0.5 * length * math.sin(math.radians(ang_deg))
    return (xc - dx, yc - dy), (xc + dx, yc + dy)


def write_case(work, name, grid, flags, source, t_end, dt_out, meta):
    d = os.path.join(work, name)
    os.makedirs(d, exist_ok=True)
    x0, y0, n, dx = grid
    with open(os.path.join(d, 'IMEX_SfloW2D.inp'), 'w') as f:
        f.write(TEMPLATE.format(name=name, t_end=t_end, dt_out=dt_out, x0=x0, y0=y0,
                                n=n, dx=dx, flags=flags, bcs=bcs(), source=source))
    flat_dem(os.path.join(d, 'flat_dem.asc'), x0, y0, n, dx)
    meta.update(dict(name=name, dx=dx, t_end=t_end, dt_out=dt_out))
    with open(os.path.join(d, 'sources.json'), 'w') as f:
        json.dump(meta, f, indent=1)


def main():
    work = sys.argv[1]
    grid = (0.0, 0.0, 100, 100.0)

    # Two fissures at different temperatures: linear ramps, a phase, a
    # second window, and one fissure oblique to the grid.
    f1 = fissure(5000.0, 5000.0, 1000.0, 0.0)
    f2 = fissure(5000.0, 3800.0, 800.0, -45.0)
    fis = [dict(p=f1, width=100.0, vel=1.0, T=800.0, tp=(30.0, 12.0, 3.0, 0.0)),
           dict(p=f2, width=100.0, vel=0.5, T=700.0, tp=(100.0, 10.0, 0.0, 95.0))]
    src = '&FISSURAL_SOURCE_PARAMETERS\n'
    for i, q in enumerate(fis, 1):
        (xa, ya), (xb, yb) = q['p']
        src += (' X_FISSURES_END_POINTS(1,%d) = %.10fD0,\n Y_FISSURES_END_POINTS(1,%d) = %.10fD0,\n'
                ' X_FISSURES_END_POINTS(2,%d) = %.10fD0,\n Y_FISSURES_END_POINTS(2,%d) = %.10fD0,\n'
                % (i, xa, i, ya, i, xb, i, yb))
        src += (' WIDTH_FISSURES(%d) = %.1fD0,\n LINEAR_VEL_FISSURES(%d) = %.1fD0,\n'
                ' T_FISSURES(%d) = %.1fD0,\n' % (i, q['width'], i, q['vel'], i, q['T']))
        for k in range(4):
            src += ' TIME_PARAM_FISSURES(%d,%d) = %.1fD0,\n' % (k + 1, i, q['tp'][k])
    src += ' ALPHAS_SOURCE = 0.02D0,\n /\n'
    write_case(work, 'fissures', grid,
               ' BOTTOM_FISSURAL_SOURCE_FLAG = T,\n N_FISSURES = 2,\n', src, 40.0, 5.0,
               dict(kind='fissures', fissures=[dict(length=1000.0, width=100.0, vel=1.0,
                                                    tp=fis[0]['tp']),
                                               dict(length=800.0, width=100.0, vel=0.5,
                                                    tp=fis[1]['tp'])]))

    # One fissure always on, at map coordinates of order 1e6 m
    gx0, gy0 = 1968000.0, 5837000.0
    (xa, ya), (xb, yb) = fissure(gx0 + 320.0, gy0 + 320.0, 160.0, 147.0)
    src = ('&FISSURAL_SOURCE_PARAMETERS\n'
           ' X_FISSURES_END_POINTS(1,1) = %.10fD0,\n Y_FISSURES_END_POINTS(1,1) = %.10fD0,\n'
           ' X_FISSURES_END_POINTS(2,1) = %.10fD0,\n Y_FISSURES_END_POINTS(2,1) = %.10fD0,\n'
           ' WIDTH_FISSURES(1) = 16.0D0,\n LINEAR_VEL_FISSURES(1) = 2.0D0,\n'
           ' T_FISSURES(1) = 800.0D0,\n ALPHAS_SOURCE = 0.02D0,\n /\n' % (xa, ya, xb, yb))
    write_case(work, 'fissure_map_coordinates', (gx0, gy0, 80, 8.0),
               ' BOTTOM_FISSURAL_SOURCE_FLAG = T,\n N_FISSURES = 1,\n', src, 10.0, 5.0,
               dict(kind='fissures', fissures=[dict(length=math.hypot(xb - xa, yb - ya),
                                                    width=16.0, vel=2.0, tp=None)]))

    # Lateral radial source: periodic window with cosine ramps, cut off
    # inside the third window by TIME_PARAM(4)
    src = ('&RADIAL_SOURCE_PARAMETERS\n X_SOURCE = 5000.0D0,\n Y_SOURCE = 5000.0D0,\n'
           ' R_SOURCE = 500.0D0,\n MFR_SOURCE = 1.0D7,\n RI_SOURCE = 0.3D0,\n'
           ' T_SOURCE = 800.0D0,\n ALPHAS_SOURCE = 0.02D0,\n'
           ' TIME_PARAM = 12.0D0, 8.0D0, 2.0D0, 30.0D0,\n /\n')
    write_case(work, 'lateral_window', grid, ' RADIAL_SOURCE_FLAG = T,\n', src, 40.0, 5.0,
               dict(kind='lateral', mfr=1.0e7, tp=(12.0, 8.0, 2.0, 30.0)))

    # Elliptical bottom source: periodic window with linear ramps and a phase
    src = ('&RADIAL_SOURCE_PARAMETERS\n X_SOURCE = 5000.0D0,\n Y_SOURCE = 5000.0D0,\n'
           ' R_SOURCE = 300.0D0,\n R2_SOURCE = 90.0D0,\n ANGLE_SOURCE = 30.0D0,\n'
           ' VEL_SOURCE = 0.5D0,\n T_SOURCE = 800.0D0,\n ALPHAS_SOURCE = 0.02D0,\n'
           ' TIME_PARAM = 25.0D0, 15.0D0, 4.0D0, 20.0D0,\n /\n')
    write_case(work, 'bottom_window', grid, ' BOTTOM_RADIAL_SOURCE_FLAG = T,\n', src, 40.0, 5.0,
               dict(kind='bottom', area=math.pi * 300.0 * 90.0, vel=0.5,
                    tp=(25.0, 15.0, 4.0, 20.0)))

    # Circular bottom source with velocity intervals
    src = ('&RADIAL_SOURCE_PARAMETERS\n X_SOURCE = 5000.0D0,\n Y_SOURCE = 5000.0D0,\n'
           ' R_SOURCE = 300.0D0,\n VEL_SOURCE = 0.5D0,\n T_SOURCE = 800.0D0,\n'
           ' ALPHAS_SOURCE = 0.02D0,\n N_INTERVALS = 3,\n'
           ' T_INTERVALS = 0.0D0, 7.3D0, 15.9D0,\n VEL_INTERVALS = 0.5D0, 0.0D0, 0.8D0,\n /\n')
    write_case(work, 'bottom_intervals', grid, ' BOTTOM_RADIAL_SOURCE_FLAG = T,\n', src, 25.0, 5.0,
               dict(kind='intervals', area=math.pi * 300.0 * 300.0,
                    t=(0.0, 7.3, 15.9), vel=(0.5, 0.0, 0.8)))


if __name__ == '__main__':
    main()
