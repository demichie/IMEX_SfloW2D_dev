Simulation example of an avalanche of finite granular mass sliding down an inclined plane and merging continuously into a horizontal plane is presented. The initial conditions and the topography of this tutorial are similar to those in Example 4.1 from Wang et al., 2004 [1]. A hemispherical shell holding the material together is suddenly released so that the bulk material commences to slide on an inclined flat plane at 35° into a horizontal run-out plane connected by a smooth transition. 
Particle sedimentation and air entrinment are neglected. 
A trasport equation for pore pressure is solved (see [2], Eq. 12), and pore pressure affects the rheological model, whic is based on the Voellmy-Salm formulation. 

A Python script (create_example.py) is provided to create the input file for this example (required packages: numpy, basemap). 
Please provide four arguments:

1) Number of cells in the x-direction 

2) Initial volume fraction of particles

3) Flow temperature (Kelvin)

4) Logical for plot of initial solution (true or false)

Usage example of the script:

>> python create_example.py 100 0.5 300 false

Run the solver:

>> ../../bin/IMEX_SfloW2D

Native NetCDF output is enabled in the input file. The simulation writes `<run_name>.nc` directly; open it with ParaView or another NetCDF-compatible tool.

REFERENCES

[1] Wang, Y., K. Hutter and S. P. Pudasaini, The Savage-Hutter theory: A system of partial differential equations for avalanche flows of snow, debris, and mud, ZAMM · Z. Angew. Math. Mech. 84, No. 8, 507 – 527 / DOI 10.1002/zamm.200310123, 2004.

[2] Gueugneau, V., Kelfoun, K., Roche, O., & Chupin, L. Effects of pore pressure in pyroclastic flows: numerical simulation and experimental validation. Geophysical Research Letters, 44(5), 2194-2202, 2017.

