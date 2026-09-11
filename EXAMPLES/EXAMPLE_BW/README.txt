This example simulate a 1d flow  with friction, deposition  and entrainment over a gentle slope. 
The total energy conservation equation is solved in this example (ENERGY_FLAG=T).
Topography does not change with deposition (TOPO_CHANGE_FLAG=F), but it is interesting to test the effect on the flow changing the flag.

The settings are similar to those presented in Bursik & Woods 1996 [1].

A Python script is provided to create the input file for this example. 
Please provide four arguments:

1) Number of cells in the x-direction 

2) Volume fraction of particles

3) Flow temperature (Kelvin)

4) Logical for plot of initial solution (true or false)

Usage example of the script:

>> ./create_example.py 400 0.01 900 false

Run the solver (this assumes that the example is in the original folder):

>> ../../bin/IMEX_SfloW2D

Native NetCDF output is enabled in the input file. The simulation writes `<run_name>.nc` directly; open it with ParaView or another NetCDF-compatible tool.


REFERENCES

[1] Bursik, M. I. &amp; Woods, A. W.
The dynamics and thermodynamics of large ash flows
Bulletin of Volcanology, 1996, 58, 175-193 
