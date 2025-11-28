from rocketcea.cea_obj_w_units import CEA_Obj

Refinement = 200

# Output options for the Solver
Graph = True
Stl = False
Dxf = True
Temperature = False
Write = False

# Engine design choices
P_combustion = 3.3 * 10 ** 6 #Pascal

obj = CEA_Obj(oxName = "N2O", fuelName = "Ethanol", temperature_units = 'degK', pressure_units='bar')
configuration = obj.get_IvacCstrTc_ChmMwGam(Pc = P_combustion / 10 ** 5, MR = 5.13, eps = 5)

L_combustion = 100 #mm
D_combustion = 50 #mm

M_exit = 2.2
g = float(configuration[4])
T_combustion = float(configuration[2])
ISP_cea = configuration[0]

R = 8.314
Mw = 0.0265 #kg/mol
Rs = R / Mw
mdot = 0.274 #kg/s

L = 5.7 #theoretical throat radius in mm

Shorten_Percentage = 0.85 #1 - Percentage to truncate nozzle by, best to have this be less than 1
Nozzle_Efficiency = 0.985 # Estimate
Combustion_Efficiency = 0.95