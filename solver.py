# Imports

import IsentropicTools as IT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
from rich import print
from rich.text import Text
from rich.rule import Rule
import CombustionChamber as CC
import Parameters as Param
import time


def m(s):
    return Text.from_markup(s)


def divider(color):
    return [Rule(style=color)]

# Main solver file for the script. Nothing in here needs to be modified by the user
# This script contains all the sub-function definitions (such as local slopes)
# And the main solver algorithm, along with the graphical output.

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
mdot = Param.mdot
M_optimal = np.sqrt(((P_combustion / 101325) ** ((g - 1) / g) - 1) * 2 / (g - 1))
L = Param.L
L_combustion = Param.L_combustion
D_combustion = Param.D_combustion
SP = Param.Shorten_Percentage
Graph2d = Param.Graph2d
Graph3d = Param.Graph3d
Materials = Param.Materials
Material = Param.Material


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
def coords_n1(k, dv, g, grid):
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
def solver(Graph2d, Graph3d, Graph3d_Fancy):

    if Graph2d: fig, ax = plt.subplots(figsize=(12, 6))
    CUSTOM_GRAY_FIG = "#161619"
    CUSTOM_GRAY_AXES = "#2C2B30"
    

    grid.set_xy(1, 1, -L / (slopes(1, 1, dv, g)[1]), 0.0)
    progress = 0
    for k1 in range(2, int(k_max)):
        x_k1, y_k1 = coords_n1(k1, dv, g, grid)
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

        if Graph2d: ax.plot(x_line, y_line, color="#D02E2E", linestyle="--", alpha=0.5)

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
        if 0 < GridField.get_x(grid, int(k_max) - NI, NI) < wall_x[-1]:
            if Graph2d: ax.plot(x_line, y_line, color="#AA5CF8", linestyle="--", alpha=0.5)
        elif GridField.get_x(grid, int(k_max) - NI, NI) >= wall_x[-1]:
            if Graph2d: ax.plot(x_line, y_line, color="#68F100", linestyle="--", alpha=0.5)
            ks.append(int(k_max) - NI + 1)
            ns.append(NI)

    y_calc = wall_y[-1]

    wall_x, wall_y = CC.CombustionChamber(wall_x, wall_y, 9.38, 47.32)
    wall_y_mirrored = -1 * wall_y
    y_min = np.min(wall_y)

    wall_y_mirrored = -wall_y

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
    Exit_Angle = np.rad2deg(
        np.arctan2(wall_y[-1] - wall_y[-2], wall_x[-1] - wall_x[-2])
    )

    time1 = time.time()

    print(f"\nComputation Time: {time1 - time0:.2f} seconds\n")

    # Plotting
    if Graph2d:

        fig.set_facecolor(CUSTOM_GRAY_FIG)
        ax.set_facecolor(CUSTOM_GRAY_AXES)

        ax.plot(
            wall_x, wall_y, color="#1E90FF", linewidth=2, label="Nozzle Wall Contour (MoC)"
        )
        ax.plot(wall_x, wall_y_mirrored, color="#1E90FF", linewidth=2)

        ax.set_title("Nozzle Wall Contour and Characteristics", color="white", fontsize=16)
        ax.set_xlabel("x (mm)", color="white")
        ax.set_ylabel("y (mm)", color="white")
        ax.axis("equal")
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.legend(
            loc="lower right",
            frameon=True,
            facecolor="black",
            edgecolor="white",
            labelcolor="white",
        )

        ax.tick_params(axis="x", colors="white")
        ax.tick_params(axis="y", colors="white")
        for spine in ax.spines.values():
            spine.set_color("white")
        ax.legend(
            loc="lower right",
            frameon=True,
            facecolor=CUSTOM_GRAY_AXES,
            edgecolor="white",
            labelcolor="white",
        )
        plt.show()
    
    if Graph3d_Fancy:

        theta = np.linspace(0, 2 * np.pi, 100)
        theta_grid, wall_x_grid = np.meshgrid(theta, wall_x)
        _, wall_y_grid = np.meshgrid(theta, wall_y)

        X = wall_x_grid
        Y = wall_y_grid * np.cos(theta_grid)
        Z = wall_y_grid * np.sin(theta_grid)

        points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

        grid2 = pv.StructuredGrid()
        grid2.points = points
        grid2.dimensions = [theta.size, wall_x.size, 1]

        grid2["Radius (mm)"] = wall_y_grid.ravel()

        plotter = pv.Plotter(window_size=[1400, 800])
        plotter.set_background("#1e1e1e") 

        plotter.add_mesh(
            grid2,
            #scalars="Radius (mm)",
            # cmap="copper",          # or 'plasma', 'turbo', 'coolwarm'
            color = "#b87333",
            opacity=0.9,
            smooth_shading=True,     # this makes a huge difference
            specular=0.4,            # metallic shine
            specular_power=120,
            ambient=0.2,
            diffuse=0.6,
            show_edges=False,
            pbr=False
        )

        profile = np.column_stack([wall_x, wall_y, np.zeros_like(wall_x)])
        profile_line = pv.Spline(profile, 300)
        plotter.add_mesh(profile_line, color="white", line_width=2)

        plotter.add_axes()
        plotter.show_grid(color="gray")

        plotter.add_title("Nozzle Geometry", font_size=14, color="white")

        plotter.show()


    if Graph3d:
        fig = plt.figure(figsize = (16,8), facecolor=CUSTOM_GRAY_FIG)
        ax = fig.add_subplot(111, projection = '3d')

        fig.set_facecolor(CUSTOM_GRAY_FIG)
        ax.set_facecolor(CUSTOM_GRAY_AXES)

        theta = np.linspace(0, 2*np.pi, 100)
        theta_grid, wall_x_grid = np.meshgrid(theta, wall_x)
        _, wall_y_grid = np.meshgrid(theta, wall_y)

        X = wall_x_grid
        Y = wall_y_grid * np.cos(theta_grid)
        Z = wall_y_grid * np.sin(theta_grid)
        norm = plt.Normalize(wall_y.min()-10, wall_y.max()+10)
        colors = cm.magma(norm(wall_y_grid))
        colors[..., 3] = 0.85

        def update_plot(high_res = True):
            ax.clear()
            ax.set_facecolor(CUSTOM_GRAY_FIG)

            # colour = Materials[Material]

            if high_res:
                ax.plot_surface(X, Y, Z, facecolors=colors, linewidth=0, antialiased=True, rstride=1, cstride=1)
            else:
                ax.plot_surface(X, Y, Z, facecolors = colors, linewidth = 0, antialiased=False, rstride = 8, cstride = 8)

        
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        ax.set_xlabel('x (mm)', color='white')
        ax.set_ylabel('y (mm)', color='white')
        ax.set_zlabel('z (mm)', color='white')
        ax.tick_params(colors='white')

        x_range = wall_x.max() - wall_x.min()
        y_range = Y.max() - Y.min()
        z_range = Z.max() - Z.min()

        plt.subplots_adjust(left=0, right=1, bottom=0, top=1)

        ax.set_xlim(wall_x.min(), wall_x.max())
        ax.set_ylim(Y.min(), Y.max())
        ax.set_zlim(Z.min(), Z.max())

        ax.set_box_aspect((x_range, y_range, z_range))

        def on_click(event):
            if event.inaxes == ax:
                update_plot(high_res=False)
        def on_release(event):
            update_plot(high_res=True)

        fig.canvas.mpl_connect('button_press_event', on_click)
        fig.canvas.mpl_connect('button_release_event', on_release)

        update_plot(high_res=True)

        plt.show()


    return (
        wall_x,
        Exit_Angle,
        y_min,
        Param.T_combustion,
        Param.g,
        P_combustion,
        M_optimal,
        IT.AreaRatio(M_optimal, g),
        wall_y,
        M_exit,
        M_exit_characteristic,
        Thrust,
        A_exit,
        P_exit,
        Ve,
        Param.ISP_cea,
        y_calc,
        L,
        wall_y_mirrored,
    )
