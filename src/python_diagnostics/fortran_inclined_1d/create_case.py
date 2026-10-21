#!/usr/bin/env python3
"""Generate the Fortran 1-D counterpart of ``kt_v2_1d.py``'s test case."""

from pathlib import Path

import numpy as np


NX = 120
DX = 0.25
X0 = 0.0
Y0 = -0.5 * DX
SLOPE = np.tan(np.deg2rad(10.0))
EXCAVATION_LEFT = 10.0
EXCAVATION_RIGHT = 20.0
EXCAVATION_DEPTH = 2.0
RHO_L = 1000.0
CP_L = 4200.0
TEMPERATURE = 300.0


def face_weight(x):
    return np.clip(
        np.minimum((x - EXCAVATION_LEFT) / DX,
                   (EXCAVATION_RIGHT - x) / DX),
        0.0, 1.0,
    )


def write_dem(path):
    # Sample the exact Python B_face function at every computational vertex.
    # Add one sample beyond both ends so the Fortran DEM interpolator has
    # support at the domain boundaries. Repeat the profile in y for 1-D flow.
    x_dem = X0 - DX + np.arange(NX + 3) * DX
    y0 = Y0
    y_dem = y0 - DX + np.arange(4) * DX
    xx, _ = np.meshgrid(x_dem, y_dem)
    z = 10.0 - SLOPE * xx - EXCAVATION_DEPTH * face_weight(xx)

    xllcorner = x_dem[0] - 0.5 * DX
    yllcorner = y_dem[0] - 0.5 * DX
    with path.open("w") as stream:
        stream.write(f"ncols {x_dem.size}\n")
        stream.write(f"nrows {y_dem.size}\n")
        stream.write(f"xllcorner {xllcorner:.12f}\n")
        stream.write(f"yllcorner {yllcorner:.12f}\n")
        stream.write(f"cellsize {DX:.12f}\n")
        stream.write("NODATA_value -9999\n")
        np.savetxt(stream, np.flipud(z), fmt="%.12f")


def write_restart(path):
    x_center = X0 + (np.arange(NX) + 0.5) * DX
    y_center = Y0 + 0.5 * DX
    x_face = X0 + np.arange(NX + 1) * DX
    weight = face_weight(x_face)
    # This is exactly the Python test's cell-average initial h.
    h = 0.5 * EXCAVATION_DEPTH * (weight[:-1] + weight[1:])
    mass = RHO_L * h
    energy = mass * CP_L * TEMPERATURE
    data = np.column_stack((x_center, np.full(NX, y_center), mass,
                            np.zeros(NX), np.zeros(NX), energy))
    with path.open("w") as stream:
        # read_solution uses fixed-width Fortran E20.12 descriptors; do not
        # put an extra delimiter between the 20-character fields.
        np.savetxt(stream, data, fmt="%20.12e", delimiter="")
        stream.write(" \n")


def write_input(path):
    path.write_text(f'''&RUN_PARAMETERS
 RUN_NAME="inclinedExcavation1D",
 RESTART=T,
 T_START=0.0D0,
 T_END=1.5D0,
 DT_OUTPUT=1.0D-1,
 OUTPUT_CONS_FLAG=T,
 OUTPUT_ESRI_FLAG=F,
 OUTPUT_NETCDF_FLAG=F,
 OUTPUT_RUNOUT_FLAG=F,
 VERBOSE_LEVEL=0,
 SERIAL_FLAG=T,
 /

&NEWRUN_PARAMETERS
 X0={X0:.12g},
 Y0={Y0:.12g},
 COMP_CELLS_X={NX},
 COMP_CELLS_Y=1,
 CELL_SIZE={DX:.12g},
 N_SOLID=0,
 N_ADD_GAS=0,
 RHEOLOGY_FLAG=F,
 GAS_FLAG=F,
 LIQUID_FLAG=T,
 RADIAL_SOURCE_FLAG=F,
 COLLAPSING_VOLUME_FLAG=F,
 TOPO_CHANGE_FLAG=F,
 SLOPE_CORRECTION_FLAG=F,
 CURVATURE_TERM_FLAG=F,
 /

&RESTART_PARAMETERS
 RESTART_FILES="inclined_excavation_1D_0000.q_2d",
 /

&WEST_BOUNDARY_CONDITIONS
 H_BCW%FLAG=1, H_BCW%VALUE=0.0D0,
 HU_BCW%FLAG=1, HU_BCW%VALUE=0.0D0,
 HV_BCW%FLAG=1, HV_BCW%VALUE=0.0D0,
 T_BCW%FLAG=1, T_BCW%VALUE=0.0D0,
 /
&EAST_BOUNDARY_CONDITIONS
 H_BCE%FLAG=1, H_BCE%VALUE=0.0D0,
 HU_BCE%FLAG=1, HU_BCE%VALUE=0.0D0,
 HV_BCE%FLAG=1, HV_BCE%VALUE=0.0D0,
 T_BCE%FLAG=1, T_BCE%VALUE=0.0D0,
 /
&SOUTH_BOUNDARY_CONDITIONS
 H_BCS%FLAG=1, H_BCS%VALUE=0.0D0,
 HU_BCS%FLAG=1, HU_BCS%VALUE=0.0D0,
 HV_BCS%FLAG=1, HV_BCS%VALUE=0.0D0,
 T_BCS%FLAG=1, T_BCS%VALUE=0.0D0,
 /
&NORTH_BOUNDARY_CONDITIONS
 H_BCN%FLAG=1, H_BCN%VALUE=0.0D0,
 HU_BCN%FLAG=1, HU_BCN%VALUE=0.0D0,
 HV_BCN%FLAG=1, HV_BCN%VALUE=0.0D0,
 T_BCN%FLAG=1, T_BCN%VALUE=0.0D0,
 /
&TEMPERATURE_PARAMETERS
 /

&NUMERIC_PARAMETERS
 KT_TRACE_RECONSTRUCTION = 'froude_blend',
 DT0=1.0D-4,
 MAX_DT=2.0D-2,
 CFL=1.4D-1,
 CFL_POLICY='directional_min',
 LIMITER=5*3,
 THETA=1.3D0,
 RECONSTR_COEFF=1.0D0,
 N_RK=2,
 /

&EXPL_TERMS_PARAMETERS
 GRAV=9.81D0,
 /

&RHEOLOGY_PARAMETERS
 RHEOLOGY_MODEL = 'none',
 /

&GAS_TRANSPORT_PARAMETERS
 SP_HEAT_A=998.0D0,
 SP_GAS_CONST_A=287.051D0,
 KIN_VISC_A=1.48D-5,
 PRES=101300.0D0,
 T_AMBIENT=300.0D0,
 ENTRAINMENT_FLAG=F,
 /

&LIQUID_TRANSPORT_PARAMETERS
 SP_HEAT_L={CP_L:.1f}D0,
 RHO_L={RHO_L:.1f}D0,
 LOSS_RATE=0.0D0,
 KIN_VISC_L=1.0D-6,
 /
''')


def main():
    case_dir = Path(__file__).resolve().parent
    write_dem(case_dir / "topography_dem.asc")
    write_restart(case_dir / "inclined_excavation_1D_0000.q_2d")
    write_input(case_dir / "IMEX_SfloW2D.inp")
    initial_face_weight = face_weight(X0 + np.arange(NX + 1) * DX)
    initial_volume = np.sum(
        0.5 * EXCAVATION_DEPTH
        * (initial_face_weight[:-1] + initial_face_weight[1:])
    ) * DX
    print(f"Generated Fortran 1-D case in {case_dir}")
    print(f"Grid: {NX} x 1, dx={DX:g} m; t_end=1.5 s, KT, n_RK=2")
    print(f"Initial volume per unit width: {initial_volume:.12g} m^2")


if __name__ == "__main__":
    main()
