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

A rectangular interior is excavated by 2 m. Each side is connected to the
surrounding plane by a continuous one-cell Q1 ramp:

    x in [10,20) m
    y in [-5,5) m.

On the interior plateau the initial water thickness is exactly 2 m. On the
four side ramps and bilinear corner patches it is h=B0-B; outside it is zero.
Therefore throughout the wet region

    eta = B + h = B0,

so the initial free surface is a plane with the same constant slope as the
original, unexcavated topography.

The slope descends toward +x (east). The physically expected initial
acceleration is therefore toward +x. There is no free-surface gradient in y.

Why the DEM format differs from old example generators
------------------------------------------------------
The padded DEM is aligned with the computational VERTICES. Its interior
samples are exactly the authoritative B_vertex values used by the continuous
HP geometry. Cell centers and shared face elevations are derived from that
single Q1 field; no discontinuous pair of bed values exists at an internal
face.

Model configuration
-------------------
  N_SOLID = 0
  N_ADD_GAS = 0
  LIQUID_FLAG = T
  GAS_FLAG = F
  RHEOLOGY_FLAG = F
  fourth equation = thermal energy
  ENTRAINMENT_FLAG = F
  LOSS_RATE = 0
  SOLVER_SCHEME = KT       (legacy input label for the sole HP-PCCU operator)
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

Expected behavior
-----------------
At t=0, u=v=0. The most important signatures are:

  1. WEST / upslope wall (x=10 m):
       water appearing for x < 10 m is spurious uphill leakage.

  2. NORTH/SOUTH walls (y=+/-5 m):
       water appearing laterally outside the excavation is also spurious,
       because the initial free-surface gradient has no y component.

  3. EAST / downslope wall (x=20 m):
       water should first accelerate inside the excavation, accumulate near
       the east ramp, and only then cross the outer crest.

The HP reconstruction uses the unique shared Q1 face bed and the
positivity-limited eta-B candidate. At each outer crest the reconstructed wet
and dry endpoint thicknesses are initially zero, so the PCCU mass transport
through all four outer faces is zero without a special step classifier.

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
