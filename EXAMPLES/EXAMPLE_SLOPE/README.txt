This example simulate a supercritical flow (Ri<1) entering from the left of the domain. An initial slope is followed by a flat topography and a discontinuity. There is no sedimentation and no entrainment.
No friction is considered in this test (RHEOLOGY_FLAG=F).
Instead of the total energy equation, a simpler transport equation (pure advection) for the temperature is solved (ENERGY_FLAG=F).

A Python script is provided to create the input file for this example. 
Please provide four arguments:

1) Number of cells in the x-direction 

2) Volume fraction of particles

3) Flow temperature (Kelvin)

4) Logical for plot of initial solution (true or false)

Usage example of the script:

>> ./create_example.py 400 1.0 300 false

Run the solver (this assumes that the example is in the original folder):

>> ../../bin/IMEX_SfloW2D

Native NetCDF output is enabled in the input file. The simulation writes `<run_name>.nc` directly; open it with ParaView or another NetCDF-compatible tool.
