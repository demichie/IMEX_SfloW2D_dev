IMEX_SfloW2D_dev - 2-D inclined excavation / wet-dry step test
================================================================

Purpose
-------
This test is designed to expose spurious numerical transport across wet/dry
interfaces on a sloping topography and to compare the current KT scheme with a
hydrostatic/well-balanced interface reconstruction.

It is based on:
  * EXAMPLES/EXAMPLE_BUMP/SUBCRITICAL for the pure-water, no-particle setup;
  * EXAMPLES/EXAMPLE_2D for the 2-D restart-file and four-boundary layout.

Geometry
--------
The unexcavated bed is a plane descending in +x:

    B0(x) = 10 - tan(10 deg) x.

Default domain:
    x = [0,30] m
    y = [-10,10] m

Default grid:
    nx = 120
    dx = dy = 0.25 m
    ny = 80

A rectangular region is excavated by a constant 2 m:

    x in [10,20) m
    y in [-5,5) m.

Inside that region the initial water thickness is exactly 2 m; outside it is
zero. Therefore inside the wet region

    eta = B + h = B0,

so the initial free surface is a plane with the same constant slope as the
original, unexcavated topography.

The slope descends toward +x (east). The physically expected initial
acceleration is therefore toward +x. There is no free-surface gradient in y.

Why the DEM format differs slightly from old example generators
---------------------------------------------------------------
The DEM is written at the computational CELL CENTERS with

    xllcorner = x_min
    yllcorner = y_min
    ncols = nx_cells
    nrows = ny_cells.

This matches the current reader coordinates
x = xllcorner + (j-0.5)*cellsize and makes the 2 m excavation jump fall exactly
between computational cells. This is useful here because the test is meant to
study a true numerical topographic step.

Model configuration
-------------------
  N_SOLID = 0
  N_ADD_GAS = 0
  LIQUID_FLAG = T
  GAS_FLAG = F
  RHEOLOGY_FLAG = F
  ENERGY_FLAG = F
  ENTRAINMENT_FLAG = F
  LOSS_RATE = 0
  SOLVER_SCHEME = KT
  CFL = 0.24
  LIMITER = generalized minmod (3)
  THETA = 1.3
  RECONSTR_COEFF = 1
  N_RK = 2
  SLOPE_CORRECTION_FLAG = F
  CURVATURE_TERM_FLAG = F

The slope/curvature corrections are deliberately disabled in this first test so
that the result can be compared directly with the simplest mathematical
analysis. A second version with SLOPE_CORRECTION_FLAG=T can be added later.

Create the default test
-----------------------
From this directory:

    python3 create_example.py 120 --plot

This creates:

    topography_dem.asc
    inclined_excavation_2D_0000.q_2d
    IMEX_SfloW2D.inp
    initial_condition_reference.npz
    initial_thickness.png          (with --plot)
    initial_centerline.png         (with --plot)

Expected behavior: current scheme
---------------------------------
At t=0, u=v=0. At the wet/dry excavation walls the dissipative part of the KT
mass flux can move water into dry cells even before the physical downslope
velocity develops.

The most important signatures are:

  1. WEST / upslope wall (x=10 m):
       water appearing for x < 10 m is spurious uphill leakage.

  2. NORTH/SOUTH walls (y=+/-5 m):
       water appearing laterally outside the excavation is also spurious,
       because the initial free-surface gradient has no y component.

  3. EAST / downslope wall (x=20 m):
       the current scheme may also leak immediately by numerical diffusion.
       Physically, with a vertical 2 m wall and zero initial velocity, the water
       should first accelerate inside the excavation, accumulate at the east
       wall, and only then overtop it.

Expected behavior: proposed hydrostatic reconstruction
------------------------------------------------------
With

    B* = max(B_L,B_R),
    h*_L = max(0, eta_L-B*),
    h*_R = max(0, eta_R-B*),

the corrected initial water depth at all four vertical excavation walls is
zero because the initial free surface equals the surrounding unexcavated bed.
Thus all four initial mass fluxes through the vertical step should be zero.

Gravity still accelerates water internally toward +x. The west/north/south
walls should remain dry. The east wall should start discharging only after the
local free surface rises above the crest.

Suggested diagnostics
---------------------
For quantitative comparison, monitor the water mass outside the original
excavation in four regions:

    M_uphill  : x < 10 m
    M_downhill: x >= 20 m
    M_south   : y < -5 m with 10 <= x < 20
    M_north   : y >= 5 m with 10 <= x < 20

For the hydrostatic scheme, M_uphill, M_south, and M_north should remain zero
(up to roundoff) until a physically driven flow actually reaches/overtops those
walls. M_downhill should remain zero at the first instant and then increase as
the downslope-moving water overtops the east wall.
