import utils.Output as Output
import engine.IsentropicTools as IT
import engine.Parameters as P
import numpy as np
import engine.solver_iteration as SI
from scipy.optimize import fsolve

##################################################################################
#                                                                                #
# Numerical Solver for Rocket Nozzles based off of the Method of Characteristics #
#                                                                                #
# Creator:  Flavio Cicero                                                        #
# Date:     28/02/2026                                                           #
# Version:  1.1                                                                  #
#                                                                                #
##################################################################################


def run(gui_mode=False):
    def engine_residuals(variables):
        mdot_guess, mach_guess = variables
        
        # Ensure we use floats to prevent indexing errors
        mdot_val = float(mdot_guess)
        mach_val = float(mach_guess)
        
        A_t = (mdot_val * P.Cstr * P.Combustion_Efficiency) / P.P_combustion
        r_t = np.sqrt(A_t / np.pi) * 1000 
        
        # SI.main returns (Cf, thrust_new, mdot_out, P_exit)
        _, thrust_new, _, P_exit = SI.main(mdot_val, r_t, mach_val)

        thrust_error = thrust_new - P.Thrust
        pressure_error = P_exit - P.Ambient_P # Perfectly expanded target

        return [thrust_error, pressure_error]

    # Preliminary estimate for initial mdot guess
    CF_in = IT.estimate_CF(P.g, 5.5, P.P_combustion)    
    mdot_initial = (P.P_combustion * (P.Thrust / (CF_in * P.P_combustion))) / P.Cstr

    initial_guess = [mdot_initial, P.M_exit]

    # Solve for perfect expansion and target thrust
    solution = fsolve(engine_residuals, x0=initial_guess)
    mdot_final, mach_final = solution

    A_t_final = (mdot_final * P.Cstr * P.Combustion_Efficiency) / P.P_combustion
    radius_throat_final = np.sqrt(A_t_final / np.pi) * 1000  

    if gui_mode:
        return {
            "mdot": mdot_final,
            "rt": radius_throat_final,
            "mach": mach_final
        }
    else:
        # Standard Terminal Output
        Output.outputTable(radius_throat_final, mdot_final, mach_final)

if __name__ == "__main__":
    run()