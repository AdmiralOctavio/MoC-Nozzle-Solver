from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.blends import newFuelBlend
from rich import print


diagnostic = False

# Mesh refinement factor. Higher = better, but slower. Default = 100.
Refinement = 100

# Output options for the Solver
Graph2d = True
Graph3d = False # Do not change
Graph3d_Fancy = False # Fancy plotting using pyvista
Stl = False
Dxf = True
Temperature = False
Write = True

# Purely for aesthetic purposes. Affects output graphs.
Materials = {
    "Copper": "#db8d5c",
    "Steel": "#525252",
    "Inconel": "#958b87",
    "Titanium": "#B4B1A7",
    "Dodger Blue": "#1E90FF",
}
Material = "Copper"

# Engine design choices
P_combustion = 3.4 * 10**6  # Pascal
Oxidiser_Fuel_Ratio = 5.13

# Fuel Definitions
ethanol80 = newFuelBlend(["Ethanol", "H2O"], [80, 20])
ethanol90 = newFuelBlend(["Ethanol", "H2O"], [90, 10])
ethanol100 = newFuelBlend(["Ethanol", "H2O"], [100, 0])

# CEA Object for specific config
obj = CEA_Obj(
    oxName="N2O", fuelName=ethanol100, temperature_units="degK", pressure_units="bar"
)
configuration = obj.get_IvacCstrTc_ChmMwGam(
    Pc=P_combustion / 10**5, MR=Oxidiser_Fuel_Ratio, eps=5
)

# Combustion Chamber dimensions
L_combustion = 93.02  # mm
Contraction_ratio = 16

# Design exit mach + other parameters
M_exit = 2.2
g = float(configuration[4])
T_combustion = float(configuration[2])
ISP_cea = configuration[0]

R = 8.314
Mw = configuration[3] / 1000  # kg/mol
Rs = R / Mw
mdot = 0.3067  # kg/s

L = 12.32 / 2  # Theoretical throat radius in mm
Chamber_Slope = 45  # Please don't change this

Shorten_Percentage = (
    0.75  # 1 - Percentage to truncate nozzle by, best to have this be less than 1
)
Nozzle_Efficiency = 0.985  # Lower bound TIC estimate.
Combustion_Efficiency = 0.85  # Also estimate.

# Diagonist to ensure CEA isn't being weird.
if diagnostic:
    Oxidisers = ["N2O", "N2O", "N2O", "N2O"]
    fuels = [ethanol80, ethanol100, ethanol80, ethanol100]
    combustion_pressures = [
        P_combustion / (10**5),
        P_combustion / (10**5),
        P_combustion / (10**5),
        P_combustion / (10**5),
    ]
    MRs = [2, 2, 5.13, 5.13]
    for i in range(len(Oxidisers)):
        obj = CEA_Obj(
            oxName=Oxidisers[i],
            fuelName=fuels[i],
            temperature_units="degK",
            pressure_units="bar",
        )
        configuration = obj.get_IvacCstrTc_ChmMwGam(
            Pc=combustion_pressures[i], MR=MRs[i], eps=5
        )

        print("[b][blue]Inputs:")
        print(
            f"[red]Oxidiser: {Oxidisers[i]}, Fuel: {fuels[i]}, Combustion Pressure: {combustion_pressures[i]:.2f} bar, Mass Ratio: {MRs[i]}"
        )
        print("[b][blue]Outputs:")
        print(
            f"Vac ISp: {configuration[0]:.2f} s, C*: {configuration[1]:.2f} m/s, Chamber Temperature: {configuration[2]:.2f} K, Gamma: {configuration[4]:.2f}"
        )
        print("\n")
