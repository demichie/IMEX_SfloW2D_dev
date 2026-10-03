#!/usr/bin/env python3
"""Generate and verify a small 2-D constant-eta lake-at-rest case."""

from pathlib import Path
import argparse
import numpy as np


NX, NY = 12, 8
DX = 1.0
X0, Y0 = 0.0, 0.0
RHO = 1000.0
CP = 4200.0
TEMP = 300.0
ETA0 = 6.0


def bed(x, y, bed_mode):
    if bed_mode == "planar":
        return 1.3 + 0.07 * x - 0.045 * y
    if bed_mode == "one-cell":
        # A steep continuous transition sampled by the input DEM.
        ramp_x = np.clip(x - 5.0, 0.0, 1.0)
        ramp_y = np.clip(y - 3.0, 0.0, 1.0)
        return 1.0 + 1.8 * ramp_x + 1.2 * ramp_y
    raise ValueError(f"unknown bed mode {bed_mode}")


def write_case(slope_correction=False, curvature=False, bed_mode="planar"):
    x = X0 + (np.arange(NX) + 0.5) * DX
    y = Y0 + (np.arange(NY) + 0.5) * DX
    X, Y = np.meshgrid(x, y)
    xd = X0 - 0.5 * DX + np.arange(NX + 2) * DX
    yd = Y0 - 0.5 * DX + np.arange(NY + 2) * DX
    Xd, Yd = np.meshgrid(xd, yd)
    dem = bed(Xd, Yd, bed_mode)

    # Match the solver's DEM -> vertex -> Q1 cell-centre construction exactly.
    bed_vertex = 0.25 * (dem[:-1, :-1] + dem[:-1, 1:] +
                         dem[1:, :-1] + dem[1:, 1:])
    bed_center = 0.25 * (bed_vertex[:-1, :-1] + bed_vertex[:-1, 1:] +
                         bed_vertex[1:, :-1] + bed_vertex[1:, 1:])
    h = ETA0 - bed_center

    initial = np.zeros((NX * NY, 6))
    initial[:, 0] = X.ravel()
    initial[:, 1] = Y.ravel()
    initial[:, 2] = (RHO * h).ravel()
    initial[:, 5] = (RHO * h * CP * TEMP).ravel()
    with Path("lake_rest_0000.q_2d").open("w") as stream:
        for row in range(NY):
            start = row * NX
            np.savetxt(stream, initial[start:start + NX], fmt="%19.12e")
            stream.write(" \n")
    np.save(f"lake_rest_reference_{bed_mode}.npy", initial)

    header = (
        f"ncols {NX + 2}\n"
        f"nrows {NY + 2}\n"
        f"xllcorner {X0 - DX}\n"
        f"yllcorner {Y0 - DX}\n"
        f"cellsize {DX}\n"
        "NODATA_value -9999\n"
    )
    np.savetxt("topography_dem.asc", np.flipud(dem), header=header,
               comments="", fmt="%.14e")

    Path("IMEX_SfloW2D.inp").write_text(f"""&RUN_PARAMETERS
 RUN_NAME="lakeRest",
 RESTART=T,
 T_START=0.0D0,
 T_END=1.0D-1,
 DT_OUTPUT=1.0D-1,
 OUTPUT_CONS_FLAG=T,
 OUTPUT_ESRI_FLAG=F,
 OUTPUT_NETCDF_FLAG=F,
 OUTPUT_RUNOUT_FLAG=F,
 VERBOSE_LEVEL=-1,
 SERIAL_FLAG=F,
 /
&NEWRUN_PARAMETERS
 X0={X0}, Y0={Y0},
 COMP_CELLS_X={NX}, COMP_CELLS_Y={NY}, CELL_SIZE={DX},
 N_SOLID=0, N_ADD_GAS=0,
 RHEOLOGY_FLAG=F, GAS_FLAG=F, LIQUID_FLAG=T,
 RADIAL_SOURCE_FLAG=F, COLLAPSING_VOLUME_FLAG=F,
 TOPO_CHANGE_FLAG=F, SLOPE_CORRECTION_FLAG={'T' if slope_correction else 'F'},
 CURVATURE_TERM_FLAG={'T' if curvature else 'F'},
 /
&RESTART_PARAMETERS
 RESTART_FILES="lake_rest_0000.q_2d",
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
 SOLVER_SCHEME="KT",
 DT0=1.0D-4, MAX_DT=1.0D-2, CFL=0.24D0,
 LIMITER=5*3, THETA=1.3D0, RECONSTR_COEFF=1.0D0,
 INTERFACES_RELAXATION=F, N_RK=2,
 /
&EXPL_TERMS_PARAMETERS
 GRAV=9.81D0,
 /
&RHEOLOGY_PARAMETERS
 RHEOLOGY_MODEL=0,
 /
&GAS_TRANSPORT_PARAMETERS
 SP_HEAT_A=998.0D0, SP_GAS_CONST_A=287.051D0,
 KIN_VISC_A=1.48D-5, PRES=101300.0D0,
 T_AMBIENT={TEMP}, ENTRAINMENT_FLAG=F,
 /
&LIQUID_TRANSPORT_PARAMETERS
 SP_HEAT_L={CP}, RHO_L={RHO}, LOSS_RATE=0.0D0,
 KIN_VISC_L=1.0D-6,
 /
""")


def check_case(output_name, bed_mode="planar"):
    reference = np.load(f"lake_rest_reference_{bed_mode}.npy")
    result = np.loadtxt(output_name)
    if result.shape != reference.shape:
        raise SystemExit(f"unexpected output shape {result.shape}")

    mass_scale = max(1.0, np.max(np.abs(reference[:, 2])))
    energy_scale = max(1.0, np.max(np.abs(reference[:, 5])))
    mass_error = np.max(np.abs(result[:, 2] - reference[:, 2])) / mass_scale
    momentum_error = np.max(np.abs(result[:, 3:5])) / mass_scale
    energy_error = np.max(np.abs(result[:, 5] - reference[:, 5])) / energy_scale
    tolerance = 2.0e-11
    if max(mass_error, momentum_error, energy_error) > tolerance:
        raise SystemExit(
            "lake-at-rest drift: "
            f"mass={mass_error:.3e}, momentum={momentum_error:.3e}, "
            f"energy={energy_error:.3e}"
        )
    print(
        f"relative errors: mass={mass_error:.3e}, "
        f"momentum={momentum_error:.3e}, energy={energy_error:.3e}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", metavar="OUTPUT")
    parser.add_argument("--slope-correction", action="store_true")
    parser.add_argument("--curvature", action="store_true")
    parser.add_argument("--bed-mode", choices=("planar", "one-cell"),
                        default="planar")
    args = parser.parse_args()
    if args.check:
        check_case(args.check, args.bed_mode)
    else:
        write_case(args.slope_correction, args.curvature, args.bed_mode)
