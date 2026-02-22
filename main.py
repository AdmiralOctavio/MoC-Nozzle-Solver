import Output
import IsentropicTools as IT
import Parameters as P
import numpy as np
import solver_iteration as SI
import scipy.optimize as opt


def run():

    # Initial thrust coefficient:
    CF_in = IT.estimate_CF(5.5, P.g, P.P_combustion)
    Cstr = P.Cstr
    A_t = P.Thrust / (CF_in * P.P_combustion)
    mdot = P.P_combustion * A_t / Cstr
    r_t = np.sqrt(A_t / np.pi) * 1000

    Cf, thrust_old, mdot, P_exit = SI.solver(mdot, r_t)

    A_t = mdot * Cstr / (P.P_combustion)
    r_t = np.sqrt(A_t / np.pi) * 1000
    M_exit_original = P.M_exit

    difference = 1

    while difference > 0.001:
        output = SI.solver(mdot, r_t)
        thrust_new = output[1]
        difference = abs((thrust_new - thrust_old) / thrust_old)
        mdot = output[2]
        A_t = mdot * Cstr / (P.P_combustion)
        r_t = np.sqrt( A_t/ np.pi) * 1000
        thrust_old = thrust_new

    radius_throat_final = np.sqrt( (output[1] / (output[0] * P.P_combustion)) / (np.pi) ) * 1000
    mdot_final = output[2]

    

    Output.outputTable(radius_throat_final, mdot_final)


if __name__ == "__main__":
    run()
