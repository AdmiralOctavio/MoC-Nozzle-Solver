from rocketcea.cea_obj_w_units import CEA_Obj

#Everything here can (and probably should) be changed to suit your own engine. Except the units. Don't touch those.

P_combustion = 3.3 * 10 ** 6 #Pascal

obj = CEA_Obj(oxName = "N2O", fuelName = "Ethanol", temperature_units = 'degK', pressure_units='bar')
configuration = obj.get_IvacCstrTc_ChmMwGam(Pc = 33, MR = 5.13, eps = 5)

L_combustion = 100 #mm
D_combustion = 50 #mm

M_exit = 2.25
g = float(configuration[4])
T_combustion = float(configuration[2])
ISP_cea = configuration[0]

R = 8.314
Mw = 0.0265 #kg/mol
Rs = R / Mw
mdot = 0.225 #kg/s

L = 5.75 #theoretical throat radius in mm

Shorten_Percentage = 0.78