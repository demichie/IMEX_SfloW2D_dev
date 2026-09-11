Simulation example of a collapse of a dilute gas-particle mixture over a topography (Campi Flegrei area). The initial conditions are similar to those presented in [1]

To run the example, unzip first the topography file.

>> unzip topo.zip

Once the topography file is unzipped, launch the solver:

>> ../../bin/IMEX_SfloW2D

Several output files are created as ESRI ascii files (*.asc) and they can be plotted with a GIS.

Native NetCDF output is enabled in the input file. The simulation writes `<run_name>.nc` directly; open it with ParaView or another NetCDF-compatible tool.


[1] A fast, calibrated model for pyroclastic density currents kinematics and hazard
TE Ongaro, S Orsucci, F Cornolti
Journal of Volcanology and Geothermal Research 327, 257-272
