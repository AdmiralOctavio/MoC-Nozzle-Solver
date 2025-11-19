import numpy as np
from scipy.optimize import brentq
from scipy.optimize import root_scalar

def PM(M, g): 
    if M >= 1:
        return np.sqrt((g+1)/(g-1)) * np.arctan(np.sqrt((g-1)/(g+1) *(M**2-1))) - np.arctan(np.sqrt(M**2-1))
    else: return 0

def mu(M):
    if M < 1.0:
        return np.nan 
    return np.arcsin(1/M) 

def PMinv(nu_target, gamma, M_guess=1.1, M_min=1.0, M_max=10.0):
    def f(M, nu_target, gamma):
        return PM(M, gamma) - nu_target
    try:
        return brentq(f, M_min, M_max, args=(nu_target, gamma))
    except ValueError:
        print(f"Warning: Could not find M for nu={nu_target}. Returning M_max.")
        return M_max 

def Pressure(P0, g, M):
    return P0 / (1 + (g-1)/2 * M**2)**(g / (g-1))

def Temperature(T0, g, M):
    return T0 / (1 + (g-1)/2 * M**2)

def LocalSoS(g, Rs, T):
    return np.sqrt(g * Rs * T)

def AreaRatio(M, g):
    term1 = (2 / (g + 1)) * (1 + (g - 1) / 2 * M**2)
    exponent = (g + 1) / (2 * (g - 1))
    return (1 / M) * term1**exponent

def AreaRatioInverse(target_AR, gamma, flow_regime):
    """
    Calculates the Mach number (M) for a given target Area Ratio (A/A*) 
    by numerically solving the transcendental equation, constrained by flow regime.
    
    Args:
        target_AR (float): The desired Area Ratio (A_exit / A_throat). Must be > 1.
        gamma (float): Specific heat ratio.
        flow_regime (str): 'subsonic' or 'supersonic' to select the root.

    Returns:
        float: The calculated Mach number M.
    """
    if target_AR < 1.0:
        raise ValueError("Area Ratio must be >= 1.0.")

    def equation_to_solve(M):
        return AreaRatio(M, gamma) - target_AR

    M_min_epsilon = 1e-6 
    M_max_supersonic = 10.0
    M_max_subsonic = 1.0 - 1e-6
    
    if flow_regime == 'supersonic':
        bracket = [1.0 + M_min_epsilon, M_max_supersonic]
    elif flow_regime == 'subsonic':
        bracket = [M_min_epsilon, M_max_subsonic]
    else:
        raise ValueError("flow_regime must be 'subsonic' or 'supersonic'.")

    try:
        solution = root_scalar(
            equation_to_solve, 
            bracket=bracket,
            method='brentq'
        )
        
        if solution.converged:
            return solution.root
        else:
            return np.nan
            
    except ValueError as e:
        print(f"Error finding root for AR={target_AR}. Check bounds. Details: {e}")
        return np.nan