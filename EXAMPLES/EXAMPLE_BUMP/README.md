# One-dimensional flow over a bump

## Purpose and role in IMEX-SfloW2D

This directory groups four independent SWASHES-inspired steady shallow-water
verification configurations. Each subdirectory contains its own input
generator, DEM and `README.md` describing the physically different
boundary conditions and expected flow regime.

## Benchmark variants

| Case | Description |
| --- | --- |
| [SUBCRITICAL](SUBCRITICAL/README.md) | Subcritical flow with prescribed inflow discharge and downstream height. |
| [SUPERCRITICAL](SUPERCRITICAL/README.md) | Supercritical flow with prescribed incoming depth/discharge and zero-gradient outflow. |
| [TRANSCRITICAL_WITHOUT_SHOCK](TRANSCRITICAL_WITHOUT_SHOCK/README.md) | Smooth subcritical-to-supercritical transition. |
| [TRANSCRITICAL_WITH_SHOCK](TRANSCRITICAL_WITH_SHOCK/README.md) | Transcritical transition including a shock. |

## Running a case

Select a variant, change into that subdirectory and follow its guide. The
`create_example.py` generators write the DEM in their local `DEM/` folder;
run the solver from that same variant's directory. Do not mix input,
topography or reference outputs across the four regimes.
