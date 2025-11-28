import IsentropicTools as IT
import numpy as np
import matplotlib.pyplot as plt
from rich import print 
import CombustionChamber as CC
import Parameters as Param
from scipy.interpolate import interp1d
from scipy.interpolate import CubicHermiteSpline
from scipy.integrate import quad
import scipy as sp
from scipy.optimize import root_scalar
import stlgenerator
import TemperatureAnalysis as TA
import time

Efficiency = Param.Nozzle_Efficiency * Param.Combustion_Efficiency
Refinement = Param.Refinement
M_exit = Param.M_exit
g = Param.g
v_e = IT.PM(M_exit, g) 
dv = (v_e - IT.PM(1.0, g)) / Refinement 
k_max = int(1/2 * v_e / dv + 1)
n_max = int(1/2 * v_e / dv + 1)
P_combustion = Param.P_combustion
T_combustion = Param.T_combustion
R = Param.R
Mw = Param.Mw
Rs = Param.Rs
mdot = Param.mdot
M_optimal = np.sqrt(((P_combustion / 101325)**((g-1)/g) - 1) * 2 / (g-1))
L = Param.L 
L_combustion = Param.L_combustion 
D_combustion = Param.D_combustion 
SP = Param.Shorten_Percentage
Graph = Param.Graph
Write = Param.Write
Stl = Param.Stl
Dxf = Param.Dxf
Temperature = Param.Temperature

fig, ax = plt.subplots(figsize=(12, 6))
CUSTOM_GRAY_FIG = "#161619"
CUSTOM_GRAY_AXES = "#2C2B30"
fig.set_facecolor(CUSTOM_GRAY_FIG) 
ax.set_facecolor(CUSTOM_GRAY_AXES) 

class GridField:
    def __init__(self, k_max, n_max):
        self.k_max = int(k_max)
        self.n_max = int(n_max)
        self.x = np.zeros((self.k_max, self.n_max))
        self.y = np.zeros((self.k_max, self.n_max))

    def set_xy(self, k, n, x_val, y_val):
        self.x[int(k)-1, int(n)-1] = x_val
        self.y[int(k)-1, int(n)-1] = y_val

    def get_xy(self, k, n):
        return self.x[int(k)-1, int(n)-1], self.y[int(k)-1, int(n)-1]

    def get_x(self, k, n):
        return self.x[int(k)-1, int(n)-1]

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

    m1 = np.tan((mu_kn + mu_km1np1)/2 + (Phi_k + Phi_km1)/2)
    m2 = -np.tan((mu_kn + mu_km1np1)/2 - (Phi_k + Phi_km1)/2)
    if k == int(k_max-n+1): return m1, np.tan(Phi_k)
    else: return m1, m2

#I think this is actually the initial kernel lowk
def coords_n1(k, dv, g, grid):
    n = 1
    x_km1, y_km1 = grid.get_xy(k - 1, n)
    m1, m2 = slopes(k, 1, dv, g)
    x_k1 = (L - (y_km1 - m1 * x_km1)) / (m1 - m2)
    y_k1 = y_km1 + m1 * (x_k1 - x_km1)
    #print(f"x_{k-1} = {x_km1}, y_{k-1} = {y_km1}, m1 = {m1}, m2 = {m2}")
    return x_k1, y_k1

def coords_k1(n, dv, g, grid):
    k = 1
    x_kp1nm1, y_kp1nm1 = grid.get_xy(k + 1, n - 1)
    m1, m2 = slopes(1, n, dv, g)
    x_kn = - (y_kp1nm1 - m2 * x_kp1nm1) / m2
    y_kn = 0
    #print(f"x_{k+1},{n-1}: {x_kp1nm1}, y_{k+1},{n-1}: {y_kp1nm1}, m2: {m2}, x_{k},{n}: {x_kn}, y_{k},{n}: {y_kn}")

    return x_kn, y_kn

def coords(k, n, dv, g, grid):
    x_kp1nm1, y_kp1nm1 = grid.get_xy(k + 1, n - 1)
    x_km1n, y_km1n = grid.get_xy(k - 1, n)

    m1, m2 = slopes(k, n, dv, g)

    x_kn = ((y_kp1nm1 - m2 * x_kp1nm1) - (y_km1n - m1 * x_km1n)) / (m1 - m2)
    y_kn = y_km1n + m1 * (x_kn - x_km1n)
    return x_kn, y_kn

time0 = time.time()

#basically the wall points are defined as kmax, 1 -> kmax - 1, 2, ... 2, nmax - 1.
def solver(Graph, Write, Model, DXF, Temperature):
    grid.set_xy(1, 1, -1/(slopes(1, 1, dv, g)[1]), 0.0)
    progress = 0
    for k1 in range(2, int(k_max) + 1):
        x_k1, y_k1 = coords_n1(k1, dv, g, grid)
        grid.set_xy(k1, 1, x_k1, y_k1)

    x_k1n2, y_k1n2 = coords_k1(2, dv, g, grid)
    grid.set_xy(1, 2, x_k1n2, y_k1n2)
    for k2 in range(2, k_max):
        x_k2, y_k2 = coords(k2, 2, dv, g, grid)
        grid.set_xy(k2, 2, x_k2, y_k2)

    for N in range(3, n_max+1):

        x_k1nN, y_k1nN = coords_k1(N, dv, g, grid)
        grid.set_xy(1, N, x_k1nN, y_k1nN)

        for kN in range(2, k_max - N +2):
            x_kN, y_kN = coords(kN, N, dv, g, grid)
            grid.set_xy(kN, N, x_kN, y_kN)
            progress += 1
        print(f" Progress: {int(progress / (n_max * k_max / 2) * 100)}% ", end='\r')

    wall_x = []
    wall_y = []

    #Prints the characteristics of Family II
    for NII in range(1, int(n_max) + 1):

        x_line = []
        y_line = []

        x_line.append(0)
        y_line.append(L)

        for NII_2 in range(1, NII+1):
            x_NN, y_NN = grid.get_xy(NII - NII_2 + 1, NII_2)
            x_line.append(x_NN)
            y_line.append(y_NN)

            if NII == n_max and NII_2 != n_max:
                wall_x.append(x_NN)
                wall_y.append(y_NN)
                
        
        ax.plot(x_line, y_line, color="#D02E2E", linestyle='--', alpha=0.5)

    split_index = int(len(wall_x) * SP)
    wall_x_truncated = wall_x[split_index:]
    wall_y_truncated = wall_y[split_index:]
    wall_x = wall_x[:split_index]
    wall_y = wall_y[:split_index]

    #Prints the characteristics of Family I
    ks = []
    ns = []
    for NI in range(1, n_max+1):
        x_line = []
        y_line = []
        for KI in range(1, int(k_max) - NI + 2):
            x_kn, y_kn = grid.get_xy(KI, NI)
            x_line.append(x_kn)
            y_line.append(y_kn)
        if 0 < GridField.get_x(grid, int(k_max) - NI + 1, NI) < wall_x[-1]:
            ax.plot(x_line, y_line, color="#AA5CF8", linestyle='--', alpha=0.5)
        elif GridField.get_x(grid, int(k_max) - NI + 1, NI) >= wall_x[-1]:  
            ax.plot(x_line, y_line, color="#68F100", linestyle='--', alpha=0.5)
            ks.append(int(k_max) - NI + 1)
            ns.append(NI)

    x_calc = grid.get_x(1, n_max-1)

    y_calc = wall_y[-1]

    Radius = wall_x[0] / np.sin(phi_k(k_max, dv))
    x_arc = np.linspace(-Radius * np.sqrt(2)/2, wall_x[0], 20)
    y_arc = -np.sqrt(Radius**2 - x_arc**2) + wall_y[0] + Radius * np.cos(phi_k(k_max, dv))

    def parabolatest():
        X = np.array([-Radius * np.sqrt(2)/2, wall_x[0]])
        Y = np.array([y_arc[0], y_arc[-1]])
        dydx = np.array([-1, np.arctan(phi_k(k_max, dv))])

        chs = CubicHermiteSpline(X, Y, dydx)

        x_poly = np.linspace(-Radius * np.sqrt(2)/2, wall_x[0], 20)
        y_poly = chs(x_poly)

        #x_poly = [0]
        #y_poly = [L]

        return x_poly, y_poly

    x_arc, y_arc = parabolatest()
    y_min = np.min(y_arc)

    wall_x = np.append(x_arc, wall_x)
    wall_y = np.append(y_arc, wall_y)

    wall_x, wall_y = CC.CombustionChamber(wall_x, wall_y, x_arc, y_arc, 30)

    wall_y_mirrored = -wall_y

    A_calc = (y_calc/1000)**2 * np.pi
    A_throat = (y_min/1000)**2 * np.pi

    x_final_characteristic, a = GridField.get_xy(grid, 1, n_max-1)
    index = np.argmin(np.abs(wall_x - x_final_characteristic))
    y_final_characteristic = wall_y[index]

    A_exit = (y_final_characteristic / 1000)**2 * np.pi
    M_exit_characteristic = IT.AreaRatioInverse(A_exit / A_throat, g, 'supersonic')

    M_exit_true = IT.AreaRatioInverse(A_calc / A_throat, g, 'supersonic')
    P_exit = IT.Pressure(P_combustion, g, M_exit_characteristic)
    T_exit = IT.Temperature(T_combustion, g, M_exit_characteristic)
    A = IT.LocalSoS(g, Rs, T_exit)

    Ve = A * M_exit_characteristic
    Thrust = (mdot * Ve + (P_exit - 101325) * A_exit) * Efficiency
    Exit_Angle = np.rad2deg(np.arctan2(wall_y[-1] - wall_y[-2], wall_x[-1] - wall_x[-2]))

    

    Temperature_profile = TA.LocalTemperature(wall_x, wall_y)


    def ExitMachDistribution(wall_x, wall_y):
        M_list = []
        x_list = []
        K_exit = np.max(ks)
        N_exit = n_max + 1 - K_exit
        CenterLine_x = np.array([])
        CenterLine_y = np.zeros(K_exit-1)
        Exit_x = np.ones(K_exit-1) * wall_x[-1]
        Exit_y = np.array([])
        Contour_x = np.array([])
        Contour_y = np.array([])
        ys = []
        for n in range(N_exit, n_max):
            x, a = GridField.get_xy(grid, 1, n)
            index = np.argmin(np.abs(wall_x - x))
            y = wall_y[index]
            ys.append(y)
            if ys[n - N_exit] == ys[n-1 - N_exit] and n > N_exit:
                ys[n - N_exit] = (ys[n - N_exit] + wall_y[index + 1])/2
            y = ys[n - N_exit]

            k = int(k_max - n + 1)
            wallx, wally = GridField.get_xy(grid, k, n)
            Contour_x = np.append(Contour_x, wallx)
            Contour_y = np.append(Contour_y, wally)
            
            CenterLine_x = np.append(CenterLine_x, x)
            M_list.append(IT.AreaRatioInverse((y**2 / y_min**2), g, 'supersonic'))
        
        def line(x1, x2, y1, y2, x):
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        
        for i in range(len(Contour_x)):
            y = line(CenterLine_x[i], Contour_x[i], 0, Contour_y[i], wall_x[-1])
            Exit_y = np.append(Exit_y, y)
        
        return np.array(M_list), Exit_y

    M_distribution, Exit_y = ExitMachDistribution(wall_x, wall_y)
    M_distribution = np.append(M_distribution, M_exit_true)
    Exit_y = np.append(Exit_y, 0)

    def MassWeightedThrust():
        coeffs = np.polyfit(Exit_y, M_distribution, 2)
        a, b, c = coeffs
        def M(y): return a * y**2 + b*y + c

        def Ve_local(y):
            return M(y) * IT.LocalSoS(g, Rs, IT.Temperature(T_combustion, g, M(y)))
        
        def Density_local(y):
            P_local = IT.Pressure(P_combustion, g, M(y))
            T_local = IT.Temperature(T_combustion, g, M(y))
            return IT.Density(P_local, Rs, T_local)

        def integrand_num(y):
            return Ve_local(y) * y * Density_local(y)
        
        def integrand_den(y):
            return y * Density_local(y)
        
        numerator, _ = quad(integrand_num, 0, wall_y[-1])
        denominator, _ = quad(integrand_den, 0, wall_y[-1])

        Ve_avg = numerator / denominator
        Thrust_mw = mdot * Ve_avg + (P_exit - 101325) * A_exit
        return Ve_avg, Thrust_mw
        
    Ve_mw, Thrust_mw = MassWeightedThrust()
    time1 = time.time()

    print(f"\nComputation Time: {time1 - time0:.2f} seconds\n")

  # --- OUTPUT SECTION ---

    print("----------------------------------------------------")
    print(f"[bold][red]Output nozzle design specifications:[/bold][/red]")
    print("[cyan3]___________________________________________________________________________________________")
    print(f"[dark_turquoise]Nozzle length: \t \t \t {wall_x[-1]:.2f} mm \t | Total length: \t \t {wall_x[-1] - wall_x[0]:.2f} mm")
    if Exit_Angle > 6: print(f"[red]Exit Angle: \t \t \t {Exit_Angle:.2f} Degrees[/red] \t [cyan3]| Exit radius: \t \t {wall_y[-1]:.2f} mm")
    else: print(f"[cyan3]Exit Angle: \t \t \t {Exit_Angle:.2f} Degrees \t | Exit radius: \t \t {wall_y[-1]:.2f} mm")
    print(f"[dark_turquoise]True Throat Radius: \t \t {y_min:.2f} mm \t | True Throat Diameter: \t {2 * y_min:.2f} mm")

    print("[dark_orange]___________________________________________________________________________________________")

    print(f"[orange_red1]Combustion Temperature:\t \t {Param.T_combustion:.2f} K \t | Gamma: \t \t \t {Param.g:.2f}")

    print("[light_sky_blue3]___________________________________________________________________________________________")

    print(f"[sky_blue2]Optimal pressure ratio: \t {P_combustion / 101325:.2f} \t \t |")
    print(f"[light_sky_blue3]Optimal exit Mach: \t \t {M_optimal:.2f} \t \t |")
    print(f"[sky_blue2]Optimal expansion ratio: \t {IT.AreaRatio(M_optimal, g):.2f} \t \t |")

    print("[light_green]___________________________________________________________________________________________")

    print(f"[green_yellow]Theoretical expansion ratio: \t {(wall_y[-1]**2 / L**2):.2f} \t \t | True expansion ratio: \t {(y_calc**2 / y_min**2):.2f}")
    print(f"[light_green]Design Exit Mach: \t \t {M_exit} \t \t | Predicted Exit Mach: \t {M_exit_characteristic:.2f}")
    print(f"[green_yellow]Predicted Thrust: \t \t {Thrust:.0f} N \t \t |")
    print(f"[light_green]> Thrust from Massflow: \t {Thrust - (P_exit - 101325) * A_exit:.2f} N \t | > Thrust from Pressure: \t {(P_exit - 101325) * A_exit:.2f} N")

    if P_exit < 0.25 * 101325: 
        print(f"[bold][red]Predicted Exit pressure: \t {P_exit:.0f} Pa")
        print(f"[bold][red]\nWARNING: Flow separation will be present at the nozzle exit.")
    elif 0.25 * 101325 < P_exit < 0.4 *101325:
        print(f"[bold][light_salmon3]Predicted Exit pressure: \t {P_exit:.0f} Pa[/light_salmon3][/bold]")
        print(f"[bold][light_salmon3]\nWARNING: Flow separation may occur at the nozzle exit.[/light_salmon3][/bold]")
    elif 0.4 * 101325 < P_exit < 0.5 * 101325: 
        print(f"[yellow3][bold]Predicted Exit pressure: \t {P_exit:.0f} Pa")
        print(f"[yellow3][bold]\nWarning: Nearing exit instability region.")
    else: print(f"[green_yellow]Predicted Exit pressure: \t {P_exit:.0f} Pa \t |")
    print(f"[light_green]Specific Impulse: \t \t {Ve / 9.80665:.2f} s \t | [light_green]Specific Impulse (CEA): \t {Param.ISP_cea:.2f} s")
    print("[blue_violet]___________________________________________________________________________________________")
    print(f"[blue_violet]Mass-Weighted Thrust: \t \t {Thrust_mw * Efficiency:.2f} N \t | [blue_violet]Mass-Weighted Isp: \t \t {Ve_mw / 9.80665 * Efficiency:.2f} s")

    ax.plot(wall_x, wall_y, color = '#1E90FF', linewidth=2, label='Nozzle Wall Contour (MoC)')
    ax.plot(wall_x, wall_y_mirrored, color = '#1E90FF', linewidth=2)

    ax.set_title("Nozzle Wall Contour and Characteristics", color='white', fontsize=16)
    ax.set_xlabel("x (mm)", color='white')
    ax.set_ylabel("y (mm)", color='white')
    ax.axis('equal')
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.legend(loc='lower right', frameon=True, facecolor='black', edgecolor='white', labelcolor='white')

    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    for spine in ax.spines.values(): spine.set_color('white')
    ax.legend(loc='lower right', frameon=True, facecolor=CUSTOM_GRAY_AXES, edgecolor='white', labelcolor='white')
    if Temperature == True:
        f_y = interp1d(wall_x, wall_y, kind='linear', bounds_error=False, fill_value="extrapolate")
        f_T = interp1d(wall_x, Temperature_profile, kind='linear', bounds_error=False, fill_value="extrapolate")

        x_plot = np.linspace(wall_x.min(), wall_x.max(), 500)
        y_plot_upper = f_y(x_plot)
        T_plot = f_T(x_plot)

        num_y_steps = 250 
        
        X_mesh, Y_mesh = np.meshgrid(x_plot, np.linspace(-1, 1, num_y_steps))

        Y_mesh_scaled = Y_mesh * y_plot_upper
        
        Z_temp = np.tile(T_plot, (num_y_steps, 1))
        
        T_min = np.min(Temperature_profile) 
        T_max = Param.T_combustion 

        gradient_fill = ax.pcolormesh(X_mesh, Y_mesh_scaled, Z_temp, 
                                    cmap='plasma', shading='auto', 
                                    vmin=T_min, vmax=T_max + 200)
        
        cbar = fig.colorbar(gradient_fill, ax=ax, orientation='vertical')
        cbar.set_label('Local Static Temperature (K)', color='white')
        cbar.ax.yaxis.set_tick_params(color='white')
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

    if Graph == True: plt.show()

    if Write == True:
        combined = np.stack((wall_x, wall_y), axis=1)
        filename = f"Nozzle_Contour_M={M_exit_true:.2f}.csv"
        np.savetxt(filename, combined, delimiter=",", header = "x_position (mm), y_radius (mm)", comments = "")

    if Model == True:
        stlgenerator.create_stl(wall_x, wall_y, M_exit_true)
    
    if DXF == True:
        stlgenerator.create_dxf(wall_x, wall_y, M_exit_true)
    
    
solver(Graph, Write, Stl, Dxf, Temperature)