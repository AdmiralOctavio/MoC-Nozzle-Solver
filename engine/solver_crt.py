# Imports

import engine.IsentropicTools as IT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
from rich import print
from rich.text import Text
from rich.rule import Rule
import engine.CombustionChamber as CC
import engine.Parameters as Param
import time

def main(mdot, L, mach):
    def m(s):
        return Text.from_markup(s)


    def divider(color):
        return [Rule(style=color)]

    # Main solver file for the script. Nothing in here needs to be modified by the user
    # This script contains all the sub-function definitions (such as local slopes)
    # And the main solver algorithm, along with the graphical output.

    # CRT Retro Terminal Palette
    CRT_BG        = "#0f0800"    # near-black warm background
    CRT_AXES      = "#100800"    # slightly lighter for axes
    CRT_AMBER     = "#ff6a00"    # primary amber
    CRT_AMBER_DIM = "#7a3000"    # dimmed amber for grids / secondary lines
    CRT_GLOW      = "#ffaa44"    # highlight / glow colour
    CRT_RED       = "#ff2200"    # char family II
    CRT_PURPLE    = "#3f36f9"    # char family I (inner)
    CRT_GREEN     = "#00ffae"    # char family I (outer)
    #CRT_BLUE      = "#31cbb4"    # wall contour
    CRT_BLUE      = "#ff6a00"


    plt.rcParams.update({
        "font.family":      "monospace",
        "font.size":        9,
        "text.color":       CRT_AMBER,
        "axes.facecolor":   CRT_AXES,
        "figure.facecolor": CRT_BG,
        "axes.edgecolor":   CRT_AMBER_DIM,
        "axes.labelcolor":  CRT_AMBER,
        "xtick.color":      CRT_AMBER,
        "ytick.color":      CRT_AMBER,
        "grid.color":       CRT_AMBER_DIM,
        "grid.linestyle":   "-",
        "grid.linewidth":   0.4,
        "grid.alpha":       0.35,
    })

    # Parameter Definition (Do not touch)

    Efficiency = Param.Nozzle_Efficiency * Param.Combustion_Efficiency
    Refinement = Param.Refinement
    M_exit = mach
    g = Param.g
    v_e = IT.PM(M_exit, g)
    dv = (v_e - IT.PM(1.0, g)) / Refinement
    k_max = int(1 / 2 * v_e / dv + 1)
    n_max = int(1 / 2 * v_e / dv + 1)
    P_combustion = Param.P_combustion
    T_combustion = Param.T_combustion
    Rs = Param.Rs
    M_optimal = np.sqrt(((P_combustion / 101325) ** ((g - 1) / g) - 1) * 2 / (g - 1))
    SP = Param.Shorten_Percentage
    Graph2d = Param.Graph2d
    Graph3d = Param.Graph3d
    R1 = Param.R1
    R2 = Param.R2
    P_amb = Param.Ambient_P


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
        if k == int(k_max - n):
            return m1, np.tan(Phi_k)
        else:
            return m1, m2


    def coords_n1(k, dv, g, grid, L):
        n = 1
        x_km1, y_km1 = grid.get_xy(k - 1, n)
        m1, m2 = slopes(k, 1, dv, g)
        x_k1 = (L - (y_km1 - m1 * x_km1)) / (m1 - m2)
        y_k1 = y_km1 + m1 * (x_k1 - x_km1)
        return x_k1, y_k1


    def coords_k1(n, dv, g, grid):
        k = 1
        x_kp1nm1, y_kp1nm1 = grid.get_xy(k + 1, n - 1)
        m1, m2 = slopes(1, n, dv, g)
        x_kn = -(y_kp1nm1 - m2 * x_kp1nm1) / m2
        y_kn = 0
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


    def solver():

        if Graph2d:
            fig, ax = plt.subplots(figsize=(14, 7))
            fig.set_facecolor(CRT_BG)
            ax.set_facecolor(CRT_AXES)

        grid.set_xy(1, 1, -L / (slopes(1, 1, dv, g)[1]), 0.0)
        progress = 0
        for k1 in range(2, int(k_max)):
            x_k1, y_k1 = coords_n1(k1, dv, g, grid, L)
            grid.set_xy(k1, 1, x_k1, y_k1)
        for N in range(2, n_max):
            x_k1nN, y_k1nN = coords_k1(N, dv, g, grid)
            grid.set_xy(1, N, x_k1nN, y_k1nN)

            for kN in range(2, k_max - N + 2):
                x_kN, y_kN = coords(kN, N, dv, g, grid)
                grid.set_xy(kN, N, x_kN, y_kN)
                progress += 1

            print(f" Progress: {int(progress / (n_max * k_max / 2) * 100)}% ", end="\r")

        wall_x = []
        wall_y = []

        # ── Helper: draw a glowing line ───────────────────────────────────────
        def glow_plot(axes, x, y, color, lw=1.5, alpha_base=1.0, zorder=3):
            """Draw a line with a soft CRT phosphor glow by layering translucent copies."""
            for w, a in [(8, 0.04), (5, 0.09), (2.5, 0.2)]:
                axes.plot(x, y, color=color, lw=lw * w / 1.5,
                          alpha=alpha_base * a, zorder=zorder - 1, solid_capstyle="round")
            axes.plot(x, y, color=color, lw=lw, alpha=alpha_base,
                      zorder=zorder, solid_capstyle="round")
            
        def glow_plot_dash(axes, x, y, color, lw=1.5, alpha_base=1.0, zorder=3, linestyle="-"):
            for w, a in [(8, 0.04), (5, 0.09), (2.5, 0.2)]:
                axes.plot(x, y, color=color, lw=lw * w / 1.5,
                        alpha=alpha_base * a, zorder=zorder - 1,
                        solid_capstyle="round", linestyle=linestyle)
            axes.plot(x, y, color=color, lw=lw, alpha=alpha_base,
                    zorder=zorder, solid_capstyle="round", linestyle=linestyle)

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

            if Graph2d:
                glow_plot(ax, x_line, y_line, CRT_RED, lw=0.8, alpha_base=0.55)

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
                if Graph2d:
                    glow_plot(ax, x_line, y_line, CRT_PURPLE, lw=0.8, alpha_base=0.8)
            elif GridField.get_x(grid, int(k_max) - NI, NI) >= wall_x[-1]:
                if Graph2d:
                    glow_plot(ax, x_line, y_line, CRT_GREEN, lw=0.8, alpha_base=0.55)
                ks.append(int(k_max) - NI + 1)
                ns.append(NI)

        y_calc = wall_y[-1]

        wall_x, wall_y = CC.CombustionChamber(wall_x, wall_y, R1, R2, L)
        wall_y_mirrored = -1 * wall_y
        y_min = np.min(wall_y)

        wall_y_mirrored = -wall_y

        A_calc = (y_calc / 1000) ** 2 * np.pi
        A_throat = (y_min / 1000) ** 2 * np.pi

        x_final_characteristic, a = GridField.get_xy(grid, 1, n_max - 1)
        index = np.argmin(np.abs(wall_x - x_final_characteristic))
        y_final_characteristic = wall_y[index]

        A_exit = (y_final_characteristic / 1000) ** 2 * np.pi
        M_exit_characteristic = IT.AreaRatioInverse(A_calc / A_throat, g, "supersonic")

        P_exit = IT.Pressure(P_combustion, g, M_exit_characteristic)
        T_exit = IT.Temperature(T_combustion, g, M_exit_characteristic)
        A = IT.LocalSoS(g, Rs, T_exit)

        Ve = A * M_exit_characteristic * Efficiency
        Thrust = (mdot * Ve + (P_exit - P_amb) * A_exit)
        Exit_Angle = np.rad2deg(
            np.arctan2(wall_y[-1] - wall_y[-2], wall_x[-1] - wall_x[-2])
        )

        time1 = time.time()

        print(f"\nComputation Time: {time1 - time0:.2f} seconds\n")

        # ── 2D Plot ───────────────────────────────────────────────────────────
        if Graph2d:

            glow_plot(ax, wall_x, wall_y,         CRT_BLUE, lw=2.0, alpha_base=1.0, zorder=5)
            glow_plot(ax, wall_x, wall_y_mirrored, CRT_BLUE, lw=2.0, alpha_base=1.0, zorder=5)
            x0 = np.linspace(np.min(wall_x), np.max(wall_x), 2)
            y0 = np.zeros(2)
            glow_plot_dash(ax, x0, y0, CRT_BLUE, lw = 1, alpha_base=0.8, zorder=5, linestyle='-.')

            y_lo = wall_y_mirrored.min() * 1.15
            y_hi = wall_y.max() * 1.15
            for scan_y in np.linspace(y_lo, y_hi, int((y_hi - y_lo) / 2)):
                ax.axhline(scan_y, color=CRT_BG, lw=0.55, alpha=0.25, zorder=6)

            ax.set_title(
                "NOZZLE WALL CONTOUR  +  CHARACTERISTICS   //   MOC SOLVER",
                color=CRT_AMBER, fontsize=11, pad=10,
            )
            ax.set_xlabel("X  (mm)", fontsize=9)
            ax.set_ylabel("Y  (mm)", fontsize=9)
            ax.axis("equal")
            ax.grid(True)

            for spine in ax.spines.values():
                spine.set_edgecolor(CRT_AMBER_DIM)
                spine.set_linewidth(1.0)

            # Legend
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color=CRT_BLUE,   lw=2,   label="NOZZLE WALL CONTOUR (MOC)"),
                Line2D([0], [0], color=CRT_RED,    lw=1,   label="CHAR. FAMILY II",   linestyle="--", alpha=0.8),
                Line2D([0], [0], color=CRT_PURPLE, lw=1,   label="CHAR. FAMILY I  (INTERNAL)", linestyle="--", alpha=0.8),
                Line2D([0], [0], color=CRT_GREEN,  lw=1,   label="CHAR. FAMILY I  (EXIT REGION)", linestyle="--", alpha=0.8),
            ]
            ax.legend(
                handles=legend_elements,
                loc="lower right",
                frameon=True,
                facecolor=CRT_AXES,
                edgecolor=CRT_AMBER_DIM,
                labelcolor=CRT_AMBER,
                fontsize=8,
            )

            # Corner label box (mimicking the terminal header)
            ax.text(
                0.01, 0.99,
                f"T-001V    MOC-SOLVER    EXIT-M: {M_exit_characteristic:.3f}",
                transform=ax.transAxes,
                color=CRT_AMBER, fontsize=7.5, va="top",
                bbox=dict(facecolor=CRT_BG, edgecolor=CRT_AMBER_DIM, lw=0.8, pad=4),
            )

        # ── 3D Plot ───────────────────────────────────────────────────────────
        if Graph3d:
            fig = plt.figure(figsize=(16, 8))
            fig.set_facecolor(CRT_BG)
            ax = fig.add_subplot(111, projection="3d")
            ax.set_facecolor(CRT_BG)

            theta = np.linspace(0, 2 * np.pi, 100)
            theta_grid, wall_x_grid = np.meshgrid(theta, wall_x)
            _, wall_y_grid = np.meshgrid(theta, wall_y)

            X = wall_x_grid
            Y = wall_y_grid * np.cos(theta_grid)
            Z = wall_y_grid * np.sin(theta_grid)

            # Amber/copper gradient mapped to radius
            norm = plt.Normalize(wall_y.min() - 10, wall_y.max() + 10)
            # Custom amber colormap: dark brown → bright amber → pale glow
            from matplotlib.colors import LinearSegmentedColormap
            amber_cmap = LinearSegmentedColormap.from_list(
                "crt_amber",
                ["#1a0800", "#7a3000", "#ff6a00", "#ffcc88"],
            )
            colors = amber_cmap(norm(wall_y_grid))
            colors[..., 3] = 0.90  # slight transparency

            def update_plot(high_res=True):
                ax.clear()
                ax.set_facecolor(CRT_BG)

                if high_res:
                    ax.plot_surface(
                        X, Y, Z,
                        facecolors=colors,
                        linewidth=0,
                        antialiased=True,
                        rstride=1, cstride=1,
                    )
                else:
                    ax.plot_surface(
                        X, Y, Z,
                        facecolors=colors,
                        linewidth=0,
                        antialiased=False,
                        rstride=8, cstride=8,
                    )

                # Glow profile line on the surface
                ax.plot(wall_x, np.zeros_like(wall_x), wall_y,
                        color=CRT_GLOW, lw=1.5, alpha=0.9, zorder=10)
                # Faint glow copy
                ax.plot(wall_x,np.zeros_like(wall_x), np.zeros_like(wall_x),
                        color=CRT_AMBER, lw=4, alpha=0.2, zorder=9)

                # Axis labels
                ax.set_xlabel("X  (mm)", color=CRT_AMBER, labelpad=8)
                ax.set_ylabel("Y  (mm)", color=CRT_AMBER, labelpad=8)
                ax.set_zlabel("Z  (mm)", color=CRT_AMBER, labelpad=8)
                ax.tick_params(colors=CRT_AMBER_DIM, labelsize=7)

                # Pane colours
                for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
                    pane.fill = True
                    pane.set_facecolor(CRT_AXES)
                    pane.set_edgecolor(CRT_AMBER_DIM)

                ax.set_title(
                    "NOZZLE GEOMETRY  //  3D SURFACE  //  MOC SOLVER",
                    color=CRT_AMBER, fontsize=11, pad=12,
                )

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

            fig.canvas.mpl_connect("button_press_event",   on_click)
            fig.canvas.mpl_connect("button_release_event", on_release)

            update_plot(high_res=True)

        print("\n" + "="*60)
        print("FINAL CALCULATION CHECK | SOLVER.PY")
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
        print(f"Pressure thrust: {(P_exit - P_amb) * A_exit:.2f} N")
        print(f"Thrust (before eff): {(mdot * Ve + (P_exit - P_amb) * A_exit):.2f} N")
        print(f"Efficiency: {Efficiency:.4f}")
        print(f"Thrust (after eff): {Thrust:.2f} N")
        print("="*60 + "\n")

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
            mdot,
            fig
        )
    return solver()