# Two-layer initial conditions

The example **no longer needs separately generated .q_2d files**.
`IMEX_SfloW2D.inp` uses the L1 thickness raster
`DEM/init_volume.asc` and the L2 dry marker `NONE`, with independent
`INITIAL_L1_PARAMETERS` and `INITIAL_L2_PARAMETERS`.

This directory is reserved for optional user-supplied conservative
snapshots, if you deliberately choose the .q_2d restart path.
It is not used by the standard Etna two-layer example.
