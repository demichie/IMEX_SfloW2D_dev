This example simulate a transcritical flow (from Ri>1, subcritical, to Ri<1, supercritical) with a shock (see [1], test 3.1.5). 
There is only one phase with constant density (1.D3 kg/s) and no gas. Entrainment, erosion and deposition are neglected in this example. 
No friction is considered in this test (RHEOLOGY_FLAG=F).
On the left boundary (subcritical) only the volumetric flow is fixed. 
On the right boundary (supercritical, right-going) zero gradient is fixed for all variables.

A Python script is provided to create the input file for this example. 
Please provide two arguments:

1) Number of cells in the x-direction 

2) Logical for plot of initial solution (true or false)

Usage example of the script:

>> ./create_example.py 200 false

Run the solver (this assumes that the example is in the original folder):

>> ../../../bin/IMEX_SfloW2D

Native NetCDF output is enabled in the input file. The simulation writes `<run_name>.nc` directly; open it with ParaView or another NetCDF-compatible tool.


REFERENCES

[1] Delestre, O., Lucas, C., Ksinant, P. A., Darboux, F., Laguerre, C., Vo, T. N. T., ... & Cordier, S. (2013). SWASHES: a compilation of shallow water analytic solutions for hydraulic and environmental studies. International Journal for Numerical Methods in Fluids, 72(3), 269-300.
