from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.blends import newFuelBlend
from rich import print
import cea
import numpy as np
from engine.Fuels import PROPELLANTS

# Default Constants & UI-Linked Parameters
P_combustion = 3.4 * 10**6
Oxidiser_Fuel_Ratio = 5.13
Refinement = 100
Thrust = 600
M_exit = 2.23
Shorten_Percentage = 0.75
Nozzle_Efficiency = 0.985
Combustion_Efficiency = 0.875
Graph2d = False
Graph3d = True
Ambient_P = 101325
Selected_Propellant = "n2o_ethanol100"  # Default

# Aesthetic / Export Options
Graph3d_Fancy = False
Stl = False
Dxf = True
Temperature = False
Write = True
Material = "Copper"
Materials = {
    "Copper": "#db8d5c",
    "Steel": "#525252",
    "Inconel": "#958b87",
    "Titanium": "#B4B1A7",
    "Dodger Blue": "#1E90FF",
}

# Physical dimensions
L_combustion = 83.02
Contraction_ratio = 16
Chamber_Slope = 45
R1 = 10
R2 = 50

# Global Containers for Solver Access
T_combustion = 0.0
Cstr = 0.0
g = 1.2
Rs = 0.0
ISP_cea = 0.0
g_exit = 1.2


def update_engine_data(pc_pa, of_ratio):
    """Calculates chemistry based on UI or manual inputs."""
    global T_combustion, Cstr, g, Rs, ISP_cea, g_exit

    # Fuel Definitions
    prop = PROPELLANTS[Selected_Propellant]
    weights = prop.get_weights(of_ratio)
    reac, prod = prop.make_mixture()
    hc = prop.get_enthalpy_coeff(of_ratio)

    solver = cea.RocketSolver(prod, reactants=reac)
    solution = cea.RocketSolution(solver)
    pc_bar = pc_pa / 1e5

    # Solve CEA
    solver.solve(
        solution, weights, pc_bar, pc_pa/Ambient_P, subar=[16], supar=[5.5], hc=hc, iac=True
    )

    T_combustion = solution.T[0]
    Cstr = solution.c_star[0]
    ISP_cea = solution.Isp[2] / 9.80665  # Using exit plane Isp

    # Get Gamma and Gas Constant using the rocketcea object for precision
    obj = prop.make_cea_obj()

    config2 = obj.get_exit_MolWt_gamma(Pc=pc_bar, MR=of_ratio, eps=5.77)
    g_exit = config2[1]
    print(f"Exit gamma: {g_exit}")

    config3 = obj.get_Throat_MolWt_gamma(Pc=pc_bar, MR=of_ratio, eps=5.77)
    g = config3[1]
    Mw = config3[0] / 1000
    Rs = 8.314 / Mw


# Initial run to populate globals
update_engine_data(P_combustion, Oxidiser_Fuel_Ratio)
