# MoC Nozzle Solver
2D Liquid Rocket Nozzle solver utilising the Method of Characteristics.

Modify relevant parameters in Parameters.py (such as desired exit mach, theoretical throat radius, etc...). Note that since this is a two dimensional solver, and not three dimensional, there will be a mismatch between the design exit mach in parameters, and the predicted exit mach. I believe this is due to the fact that the 2D solver approximates the contour as a continuous plane, and not as a cylindrical structure. Nonetheless, this tool has been verified by comparing with contours produced with RPA-4's TIC Nozzle generator, which utilises the 3D MoC method. There is a slight deviation in geometry, but generally within 1%. Why use this over RPA-4? This is free :).

Run main.py for tabular and graphical output. 

This method is based off of NASA Technical Memorandum TM X-1502 "Computer program for design of two-dimensional supersonic nozzle with sharp-edged throat", modified for use in Python along with other modifications.

Due to the iterative nature of the Method of Characteristics, and physical limitations with the design, you will always notice a discrepancy between the design exit mach, and the predicted exit mach. This is to be expected, as the two are calculated in fundamentally different ways. The design exit mach is calculated at the grid point 1, n_max - 1, whereas the predicted exit mach is calculated utilising the nozzle area ratio. This does not mean that the nozzle is non-isentropic however - the nozzle is still (theoretically) fully isentropic, but for your predicted exit Mach, not your design exit Mach. 


## This software requires the gfortran compiler for rocketCEA. There are two main ways to get this working quickly and smoothly.

### Option A (Recommended):
1 - Install Miniconda

2 - Open Anaconda Prompt, and cd into the directory where you've installed the script
```bash
cd C:\Users\XXXX\Documents\GitHub\MoC-Nozzle-Solver
```
3 - To create your virtual environment, run:
```bash
conda env create -f environment.yml
```
4 - to activate your virtual environment, run:
```bash
conda activate nozzle-env
```
5 - To install the dependencies, run:
```bash
pip install -e .
```
6 - To actually run the script, run:
```bash
run-nozzle
```

### Option B:
You can also run this using standard python (pip only), but it can be a bit more annoying.
You must manually install the fortran compiler (I think its like MinGW-w64 or sm, with gfortran selected)

1 - Set up virtual environment using 
```bash
python -m venv venv.\venv\Scripts\activate
```
2 - install the project using 
```bash
pip install -e .
```
3 - Run the solver with run-nozzle

Do this in your IDE terminal environment
