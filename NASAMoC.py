import IsentropicTools as IT
import numpy as np
import matplotlib.pyplot as plt
from rich import print 
import CombustionChamber as CC
import Parameters as Param

M_exit = Param.M_exit
g = Param.g
v_e = IT.PM(M_exit, g) #Exit Prandtl-Meyer angle
dv = (v_e - IT.PM(1.0, g)) / 20 #Incremental angle change

P_combustion = Param.P_combustion
T_combustion = Param.T_combustion

R = Param.R
Mw = Param.Mw
Rs = Param.Rs
mdot = Param.mdot

k_max = int(1/2 * v_e / dv + 1)
n_max = int(1/2 * v_e / dv + 1)

L = Param.L #theoretical throat radius in mm
L_combustion = Param.L_combustion #mm
D_combustion = Param.D_combustion #mm
SP = Param.Shorten_Percentage

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

#basically the wall points are defined as kmax, 1 -> kmax - 1, 2, ... 2, nmax - 1.

grid.set_xy(1, 1, -1/(slopes(1, 1, dv, g)[1]), 0.0)

#for K in range(1, int(k_max) + 1):
#    grid.set_xy(K, 1, 0.0, L)

print(f"\n Max N, K = {n_max}")

print("\n Initial Kernel Starting \n")
print(f"Set point k={1}, n=1 to x={-1/(slopes(1, 1, dv, g)[0])}, y=0")
for k1 in range(2, int(k_max) + 1):
    x_k1, y_k1 = coords_n1(k1, dv, g, grid)
    grid.set_xy(k1, 1, x_k1, y_k1)
    print(f"Set point k={k1}, n=1 to x={x_k1}, y={y_k1}")
print("\n Initial Kernel Complete \n")

print("\n Beginning second layer \n")

x_k1n2, y_k1n2 = coords_k1(2, dv, g, grid)
grid.set_xy(1, 2, x_k1n2, y_k1n2)
print(f"Set point k={1}, n=2 to x={x_k1n2}, y={y_k1n2}")
for k2 in range(2, k_max):
    x_k2, y_k2 = coords(k2, 2, dv, g, grid)
    grid.set_xy(k2, 2, x_k2, y_k2)
    print(f"Set point k={k2}, n=2 to x={x_k2}, y={y_k2}")

print("\n Second Layer Complete \n")

for N in range(3, n_max+1):

    print(f"\n Beginning Layer No. {N}\n")

    x_k1nN, y_k1nN = coords_k1(N, dv, g, grid)
    grid.set_xy(1, N, x_k1nN, y_k1nN)

    for kN in range(2, k_max - N +2):
        x_kN, y_kN = coords(kN, N, dv, g, grid)
        grid.set_xy(kN, N, x_kN, y_kN)
        print(f"Set point k={kN}, n={N} to x={x_kN}, y={y_kN}")
    
    print(f"\n Ending Layer No. {N}\n")


#Prints the characteristics of Family I
for NI in range(1, n_max+1):
    x_line = []
    y_line = []

    for KI in range(1, int(k_max) - NI + 2):
        x_kn, y_kn = grid.get_xy(KI, NI)
        x_line.append(x_kn)
        y_line.append(y_kn)
    plt.plot(x_line, y_line, color='darkgray', linestyle='--')

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
            
    
    plt.plot(x_line, y_line, color='lightblue', linestyle='--')

split_index = int(len(wall_x) * SP)
wall_x_truncated = wall_x[split_index:]
wall_y_truncated = wall_y[split_index:]
wall_x = wall_x[:split_index]
wall_y = wall_y[:split_index]



# Calculating the circular radius near the throat, because I don't like how NASA just has a sharp edge + circular based 


Radius = wall_x[0] / np.sin(phi_k(k_max, dv))
# at what x position would it be 45 degrees?
# R / sqrt(2) tbh
x_arc = np.linspace(-Radius * np.sqrt(2)/2, wall_x[0], 20)
y_arc = -np.sqrt(Radius**2 - x_arc**2) + wall_y[0] + Radius * np.cos(phi_k(k_max, dv))

wall_x = np.append(x_arc, wall_x)
wall_y = np.append(y_arc, wall_y)

wall_x, wall_y = CC.CombustionChamber(wall_x, wall_y, x_arc, y_arc)

wall_y_mirrored = -wall_y
#plt.plot(x_arc, y_arc, color = 'blue')

x_calc, y_calc = grid.get_xy(1, n_max-1)

difference = np.abs(wall_x - x_calc)
index = difference.argmin()
if wall_x[index] - x_calc < 0:
    y_calc = wall_y[index] + (wall_y[index+1] - wall_y[index])/(wall_x[index+1] - wall_x[index]) * (x_calc - wall_x[index])
else:
    y_calc = wall_y[index] + (wall_y[index] - wall_y[index-1])/(wall_x[index] - wall_x[index-1]) * (x_calc - wall_x[index-1])


A_calc = (y_calc/1000)**2 * np.pi
A_throat = (y_arc[5]/1000)**2 * np.pi

M_exit_true = IT.AreaRatioInverse(A_calc / A_throat, g)
P_exit = IT.Pressure(P_combustion, g, M_exit_true)
T_exit = IT.Temperature(T_combustion, g, M_exit_true)
A = IT.LocalSoS(g, Rs, T_exit)
A_exit = (wall_y[-1] / 1000)**2 * np.pi
Ve = A * M_exit_true
Thrust = mdot * Ve + (P_exit - 101325) * A_exit
Exit_Angle = np.rad2deg(np.arctan2(wall_y[-1] - wall_y[-2], wall_x[-1] - wall_x[-2]))

print(f"[bold]Output nozzle design specifications:[/bold]")
print(f"Total length: {wall_x[-1]:.2f} mm")
print(f"Exit radius: {wall_y[-1]:.2f} mm")
print(f"Exit Angle: {Exit_Angle:.2f} Degrees")
print(f"Throat Radius: {y_arc[5]:.2f} mm \n")
print(f"Theoretical expansion ratio: {(wall_y[-1]**2 / L**2):.2f}")
print(f"True expansion ratio: {(wall_y[-1]**2 / y_arc[5]**2):.2f}")
print(f"Design Exit Mach: {M_exit}")
print(f"Predicted Exit Mach: {M_exit_true:.2f}")
print(f"Predicted Thrust: {Thrust:.2f} N")
print(f"Predicted Exit pressure: {P_exit} Pa")

plt.plot(wall_x, wall_y, color = 'blue')
plt.plot(wall_x, wall_y_mirrored, color = 'blue')
plt.title("Nozzle Wall Contour")
plt.xlabel("x (mm)")
plt.ylabel("y (mm)")
plt.axis('equal')
plt.grid(True)
plt.show()