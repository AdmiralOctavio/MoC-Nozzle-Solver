import numpy as np
from scipy.optimize import brentq
from scipy.optimize import root_scalar


def PM(M, g):
    if M >= 1:
        return np.sqrt((g + 1) / (g - 1)) * np.arctan(
            np.sqrt((g - 1) / (g + 1) * (M**2 - 1))
        ) - np.arctan(np.sqrt(M**2 - 1))
    else:
        return 0


def mu(M):
    if M < 1.0:
        return np.nan
    return np.arcsin(1 / M)


def PMinv(nu_target, gamma, M_guess=1.1, M_min=1.0, M_max=10.0):
    def f(M, nu_target, gamma):
        return PM(M, gamma) - nu_target

    try:
        return brentq(f, M_min, M_max, args=(nu_target, gamma))
    except ValueError:
        print(f"Warning: Could not find M for nu={nu_target}. Returning M_max.")
        return M_max


def Pressure(P0, g, M):
    return P0 / (1 + (g - 1) / 2 * M**2) ** (g / (g - 1))


def Temperature(T0, g, M):
    return T0 / (1 + (g - 1) / 2 * M**2)


def LocalSoS(g, Rs, T):
    return np.sqrt(g * Rs * T)


def Density(P, Rs, T):
    return P / (Rs * T)


def AreaRatio(M, g):
    term1 = (2 / (g + 1)) * (1 + (g - 1) / 2 * M**2)
    exponent = (g + 1) / (2 * (g - 1))
    return (1 / M) * term1**exponent


def AreaRatioInverse(target_AR, gamma, flow_regime):
    if target_AR < 1.0:
        raise ValueError("Area Ratio must be >= 1.0.")

    def equation_to_solve(M):
        return AreaRatio(M, gamma) - target_AR

    M_min_epsilon = 1e-6
    M_max_supersonic = 10.0
    M_max_subsonic = 1.0 - 1e-6

    if flow_regime == "supersonic":
        bracket = [1.0 + M_min_epsilon, M_max_supersonic]
    elif flow_regime == "subsonic":
        bracket = [M_min_epsilon, M_max_subsonic]
    else:
        raise ValueError("flow_regime must be 'subsonic' or 'supersonic'.")

    try:
        solution = root_scalar(equation_to_solve, bracket=bracket, method="brentq")

        if solution.converged:
            return solution.root
        else:
            return np.nan

    except ValueError as e:
        print(f"Error finding root for AR={target_AR}. Check bounds. Details: {e}")
        return np.nan


def estimate_CF(gamma, epsilon, Pc):

    Me = AreaRatioInverse(epsilon, gamma, "supersonic")

    Pamb = 101325
    Pe_Pc = (1 + (gamma - 1) / 2 * Me**2) ** (-gamma / (gamma - 1))
    Pe = Pe_Pc * Pc
    
    term1 = np.sqrt((2 * gamma**2) / (gamma - 1))
    term2 = (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1))
    term3 = 1 - Pe_Pc ** ((gamma - 1) / gamma)
    CF_momentum = np.sqrt(term1 * term2 * term3)
    
    CF_pressure = epsilon * (Pe - Pamb) / Pc
    
    CF_total = CF_momentum + CF_pressure
    
    return CF_total