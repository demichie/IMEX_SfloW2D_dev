This example simulate a 1D dam-break problem with deposition, entrainment and friction over a flat topography. 

A Python script is provided to create the input file for this example. 
Please provide four arguments:

1) Number of cells in the x-direction 

2) Volume fraction of particles

3) Flow temperature (Kelvin)

4) Logical for plot of initial solution (true or false)

Usage example of the script:

>> ./create_example.py 400 0.5 300 false

Run the solver:

>> ../../bin/IMEX_SfloW2D

Native NetCDF output is enabled in the input file. The simulation writes `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc`; open it with ParaView or another NetCDF-compatible tool.
