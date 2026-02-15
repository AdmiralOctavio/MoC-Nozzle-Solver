# MoC Nozzle Solver
Liquid Rocket Nozzle solver utilising the method of characteristics.

Modify relevant parameters in Parameters.py (such as desired exit mach, theoretical throat radius, etc...).

Run main.py for tabular and graphical output. 

This method is based off of NASA Technical Memorandum TM X-1502 "Computer program for design of two-dimensional supersonic nozzle with sharp-edged throat", modified for use in Python along with other modifications.

Due to the iterative nature of the Method of Characteristics, and physical limitations with the design, you will always notice a discrepancy between the design exit mach, and the predicted exit mach. This is to be expected, as the two are calculated in fundamentally different ways. The design exit mach is calculated at the grid point 1, n_max - 1, whereas the predicted exit mach is calculated utilising the nozzle area ratio. This does not mean that the nozzle is non-isentropic however - the nozzle is still (theoretically) fully isentropic, but for your predicted exit Mach, not your design exit Mach. 

Remember to manually modify the inputs to the main solver function "solver" in NASAMoC.py to get the desired outputs (Graphs, Temperature gradient, dxf files, etc).

Dependencies:
Numpy
Scipy
Matplotlib
Rich
RocketCEA
Fortran Compiler
numpy-stl
ezdxf
os