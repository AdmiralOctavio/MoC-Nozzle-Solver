# Imports

import IsentropicTools as IT
import numpy as np
import CombustionChamber as CC
import Parameters as Param
import time

# Iterative solver file for the script. Nothing in here needs to be modified by the user
# This script contains all the sub-function definitions (such as local slopes)
# And the iterative solver algorithm.

# Parameter Definition (Do not touch)

Efficiency = Param.Nozzle_Efficiency * Param.Combustion_Efficiency
Refinement = Param.Refinement
M_exit = Param.M_exit
g = Param.g
v_e = IT.PM(M_exit, g)
dv = (v_e - IT.PM(1.0, g)) / Refinement
k_max = int(1 / 2 * v_e / dv + 1)
n_max = int(1 / 2 * v_e / dv + 1)
P_combustion = Param.P_combustion
T_combustion = Param.T_combustion
R = Param.R
Mw = Param.Mw
Rs = Param.Rs
M_optimal = np.sqrt(((P_combustion / 101325) ** ((g - 1) / g) - 1) * 2 / (g - 1))
L_combustion = Param.L_combustion
SP = Param.Shorten_Percentage
R1 = Param.R1
R2 = Param.R2
Cstr = Param.Cstr
Thrust_target = Param.Thrust

class GridField:
    def __init__(self, k_max, n_max):
        self.k_max = int(k_max)
        self.n_max = int(n_max)
        self.M = np.ones((self.k_max, self.n_max))
        self.x = np.zeros((self.k_max, self.n_max))
        self.y = np.zeros((self.k_max, self.n_max))

    def set_xy(self, k, n, x_val, y_val):
        self.x[int(k) - 1, int(n) - 1] = x_val
        self.y[int(k) - 1, int(n) - 1] = y_val

    def set_M(self, k, n, M_val):
        self.M[int(k) - 1, int(n) - 1] = M_val

    def get_M(self, k, n):
        return self.M[int(k) - 1, int(n) - 1]

    def get_xy(self, k, n):
        return self.x[int(k) - 1, int(n) - 1], self.y[int(k) - 1, int(n) - 1]

    def get_x(self, k, n):
        return self.x[int(k) - 1, int(n) - 1]

    def as_arrays(self):
        return self.x, self.y


grid = GridField(k_max, n_max)


def phi_k(k, dv):
    return (k - 1) * dv


def v_region(k, n, dv):
    return 2 * dv * (n - 1) + (k - 1) * dv


def slopes(k, n, dv, g):
    v_kn = v_region(k, n, dv)
    v_km1np1 = v_region(k - 1, n + 1, dv)

    M_kn = IT.PMinv(v_kn, g)
    M_km1np1 = IT.PMinv(v_km1np1, g)

    mu_kn = IT.mu(M_kn)
    mu_km1np1 = IT.mu(M_km1np1)

    Phi_k = phi_k(k, dv)
    Phi_km1 = phi_k(k - 1, dv)

    m1 = np.tan((mu_kn + mu_km1np1) / 2 + (Phi_k + Phi_km1) / 2)
    m2 = -np.tan((mu_kn + mu_km1np1) / 2 - (Phi_k + Phi_km1) / 2)
    # print(f"\n Wall k: {k_max - n + 1}")
    if k == int(k_max - n):
        return m1, np.tan(Phi_k)
    else:
        return m1, m2


# I think this is actually the initial kernel lowk
def coords_n1(k, dv, g, grid, L):
    n = 1
    x_km1, y_km1 = grid.get_xy(k - 1, n)
    m1, m2 = slopes(k, 1, dv, g)
    x_k1 = (L - (y_km1 - m1 * x_km1)) / (m1 - m2)
    y_k1 = y_km1 + m1 * (x_k1 - x_km1)
    # print(f"x_{k-1} = {x_km1}, y_{k-1} = {y_km1}, m1 = {m1}, m2 = {m2}")
    return x_k1, y_k1


def coords_k1(n, dv, g, grid):
    k = 1
    x_kp1nm1, y_kp1nm1 = grid.get_xy(k + 1, n - 1)
    m1, m2 = slopes(1, n, dv, g)
    x_kn = -(y_kp1nm1 - m2 * x_kp1nm1) / m2
    y_kn = 0
    # print(f"x_{k+1},{n-1}: {x_kp1nm1}, y_{k+1},{n-1}: {y_kp1nm1}, m2: {m2}, x_{k},{n}: {x_kn}, y_{k},{n}: {y_kn}")

    return x_kn, y_kn


def coords(k, n, dv, g, grid):
    x_kp1nm1, y_kp1nm1 = grid.get_xy(k + 1, n - 1)
    x_km1n, y_km1n = grid.get_xy(k - 1, n)

    m1, m2 = slopes(k, n, dv, g)

    x_kn = ((y_kp1nm1 - m2 * x_kp1nm1) - (y_km1n - m1 * x_km1n)) / (m1 - m2)
    y_kn = y_km1n + m1 * (x_kn - x_km1n)
    return x_kn, y_kn


time0 = time.time()
print(n_max * k_max / 2)


# basically the wall points are defined as kmax, 1 -> kmax - 1, 2, ... 2, nmax - 1.
def solver(mdot, L):



    grid.set_xy(1, 1, -L / (slopes(1, 1, dv, g)[1]), 0.0)
    progress = 0
    for k1 in range(2, int(k_max)):
        x_k1, y_k1 = coords_n1(k1, dv, g, grid, L)
        grid.set_xy(k1, 1, x_k1, y_k1)
        # print(f"K: {k1}, N = 1, x: {x_k1}, y: {y_k1}")
    for N in range(2, n_max):
        x_k1nN, y_k1nN = coords_k1(N, dv, g, grid)
        grid.set_xy(1, N, x_k1nN, y_k1nN)
        # print(f"K: 1, N = {N}, x: {x_k1nN}, y: {y_k1nN}")

        for kN in range(2, k_max - N + 2):
            x_kN, y_kN = coords(kN, N, dv, g, grid)
            grid.set_xy(kN, N, x_kN, y_kN)
            # print(f"K: {kN}, N = {N}, x: {x_kN}, y: {y_kN}")
            progress += 1

        print(f" Progress: {int(progress / (n_max * k_max / 2) * 100)}% ", end="\r")

    wall_x = []
    wall_y = []

    # Prints the characteristics of Family II
    for NII in range(1, int(n_max)):
        x_line = []
        y_line = []

        x_line.append(0)
        y_line.append(L)

        for NII_2 in range(1, NII + 1):
            x_NN, y_NN = grid.get_xy(NII - NII_2 + 1, NII_2)
            x_line.append(x_NN)
            y_line.append(y_NN)

            if NII == n_max - 1 and NII_2 != n_max - 1:
                wall_x.append(x_NN)
                wall_y.append(y_NN)

    split_index = int(len(wall_x) * SP)
    wall_x = wall_x[:split_index]
    wall_y = wall_y[:split_index]

    # Prints the characteristics of Family I
    ks = []
    ns = []
    for NI in range(1, n_max - 1):
        x_line = []
        y_line = []
        for KI in range(1, int(k_max) - NI + 1):
            x_kn, y_kn = grid.get_xy(KI, NI)
            x_line.append(x_kn)
            y_line.append(y_kn)
        if GridField.get_x(grid, int(k_max) - NI, NI) >= wall_x[-1]:
            ks.append(int(k_max) - NI + 1)
            ns.append(NI)

    y_calc = wall_y[-1]

    wall_x, wall_y = CC.CombustionChamber(wall_x, wall_y, R1, R2, L)

    y_min = np.min(wall_y)


    A_calc = (y_calc / 1000) ** 2 * np.pi
    A_throat = (y_min / 1000) ** 2 * np.pi

    x_final_characteristic, a = GridField.get_xy(grid, 1, n_max - 1)
    index = np.argmin(np.abs(wall_x - x_final_characteristic))
    y_final_characteristic = wall_y[index]

    A_exit = (y_final_characteristic / 1000) ** 2 * np.pi
    # M_exit_characteristic = IT.AreaRatioInverse(A_exit / A_throat, g, 'supersonic')
    M_exit_characteristic = IT.AreaRatioInverse(A_calc / A_throat, g, "supersonic")

    P_exit = IT.Pressure(P_combustion, g, M_exit_characteristic)
    T_exit = IT.Temperature(T_combustion, g, M_exit_characteristic)
    A = IT.LocalSoS(g, Rs, T_exit)

    Ve = A * M_exit_characteristic
    Thrust = (mdot * Ve + (P_exit - 101325) * A_exit) * Efficiency

    time1 = time.time()

    print(f"\nComputation Time: {time1 - time0:.2f} seconds\n")

    C_F = Thrust / (P_combustion * A_throat)
    A_throat_new = (Thrust_target) / (C_F * P_combustion)

    mdot_new = P_combustion * A_throat_new / Cstr

    print("\n" + "="*60)
    print("FINAL CALCULATION CHECK | SOLVER_ITERATION.PY")
    print("="*60)
    print(f"Input L (throat radius): {L:.4f} mm")
    print(f"Input mdot: {mdot:.6f} kg/s")
    print(f"Calculated y_min: {y_min:.4f} mm")
    print(f"A_throat: {A_throat*1e6:.4f} mm²")
    print(f"A_exit: {A_exit*1e6:.4f} mm²")
    print(f"Expansion ratio: {A_exit/A_throat:.4f}")
    print(f"M_exit: {M_exit_characteristic:.4f}")
    print(f"P_exit: {P_exit/1e5:.4f} bar")
    print(f"T_exit: {T_exit:.2f} K")
    print(f"Ve: {Ve:.2f} m/s")
    print(f"Momentum thrust: {mdot * Ve:.2f} N")
    print(f"Pressure thrust: {(P_exit - 101325) * A_exit:.2f} N")
    print(f"Thrust (before eff): {(mdot * Ve + (P_exit - 101325) * A_exit):.2f} N")
    print(f"Efficiency: {Efficiency:.4f}")
    print(f"Thrust (after eff): {Thrust:.2f} N")
    print(f"C_F: {C_F:.4f}")
    print("="*60 + "\n")

    return C_F, Thrust, mdot_new, P_exit