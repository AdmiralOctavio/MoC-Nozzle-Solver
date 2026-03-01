# Imports

import engine.IsentropicTools as IT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
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

    # ── LCARS Palette ─────────────────────────────────────────────────────────
    LCARS_BG = "#0d0a08"  # deep navy black
    LCARS_SURFACE = "#0d0a08"  # panel surface
    LCARS_ORANGE = "#ff9900"  # primary LCARS orange
    LCARS_GOLD = "#ffcc66"  # secondary gold
    LCARS_PURPLE = "#9977cc"  # LCARS purple
    LCARS_TEAL = "#44aaaa"  # LCARS teal
    LCARS_RED = "#cc4444"  # warning red
    LCARS_DIM = "#444466"  # dimmed grid/border

    # ── Global rcParams ───────────────────────────────────────────────────────
    plt.rcParams.update(
        {
            "font.family": "monospace",
            "font.size": 9,
            "text.color": LCARS_GOLD,
            "axes.facecolor": LCARS_SURFACE,
            "figure.facecolor": LCARS_BG,
            "axes.edgecolor": LCARS_DIM,
            "axes.labelcolor": LCARS_GOLD,
            "xtick.color": LCARS_GOLD,
            "ytick.color": LCARS_GOLD,
            "grid.color": LCARS_DIM,
            "grid.linestyle": "-",
            "grid.linewidth": 0.4,
            "grid.alpha": 0.4,
        }
    )

    # ── Parameter Definition ──────────────────────────────────────────────────
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
        return m1, m2

    def coords_n1(k, dv, g, grid, L):
        x_km1, y_km1 = grid.get_xy(k - 1, 1)
        m1, m2 = slopes(k, 1, dv, g)
        x_k1 = (L - (y_km1 - m1 * x_km1)) / (m1 - m2)
        y_k1 = y_km1 + m1 * (x_k1 - x_km1)
        return x_k1, y_k1

    def coords_k1(n, dv, g, grid):
        x_kp1nm1, y_kp1nm1 = grid.get_xy(2, n - 1)
        m1, m2 = slopes(1, n, dv, g)
        x_kn = -(y_kp1nm1 - m2 * x_kp1nm1) / m2
        return x_kn, 0.0

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

        # ── 2D figure setup ───────────────────────────────────────────────────
        if Graph2d:
            fig = plt.figure(figsize=(15, 7))
            fig.set_facecolor(LCARS_BG)

            # ── Left rail (colour blocks) ─────────────────────────────────────
            ax_rail = fig.add_axes([0.0, 0.0, 0.055, 1.0])
            ax_rail.set_facecolor(LCARS_BG)
            ax_rail.axis("off")

            rail_colors = [
                LCARS_ORANGE,
                LCARS_PURPLE,
                LCARS_TEAL,
                LCARS_GOLD,
                LCARS_ORANGE,
                LCARS_PURPLE,
                LCARS_TEAL,
            ]
            rail_heights = [0.16, 0.05, 0.20, 0.05, 0.24, 0.05, 0.16]
            y_cur = 0.97
            gap = 0.012
            for rc, rh in zip(rail_colors, rail_heights):
                ax_rail.add_patch(
                    plt.Rectangle(
                        (0.15, y_cur - rh),
                        0.7,
                        rh,
                        transform=ax_rail.transAxes,
                        color=rc,
                        clip_on=False,
                    )
                )
                y_cur -= rh + gap

            # Rounded top cap on rail
            from matplotlib.patches import FancyBboxPatch

            ax_rail.add_patch(
                FancyBboxPatch(
                    (0.15, 0.955),
                    0.7,
                    0.04,
                    boxstyle="round,pad=0.015",
                    transform=ax_rail.transAxes,
                    facecolor=LCARS_ORANGE,
                    edgecolor="none",
                    clip_on=False,
                )
            )

            # ── Top header bar ────────────────────────────────────────────────
            ax_top = fig.add_axes([0.07, 0.88, 0.92, 0.075])
            ax_top.set_facecolor(LCARS_BG)
            ax_top.axis("off")

            header_segs = [
                (0.000, 0.005, LCARS_BG),
                (0.005, 0.480, LCARS_ORANGE),
                (0.487, 0.004, LCARS_BG),
                (0.491, 0.130, LCARS_TEAL),
                (0.623, 0.004, LCARS_BG),
                (0.627, 0.180, LCARS_PURPLE),
                (0.809, 0.004, LCARS_BG),
                (0.813, 0.187, LCARS_GOLD),
            ]
            for sx, sw, sc in header_segs:
                ax_top.add_patch(
                    plt.Rectangle(
                        (sx, 0),
                        sw,
                        1.0,
                        transform=ax_top.transAxes,
                        color=sc,
                    )
                )
            ax_top.text(
                0.01,
                0.5,
                "NOZZLE WALL CONTOUR  +  CHARACTERISTICS   //   MOC SOLVER",
                transform=ax_top.transAxes,
                color=LCARS_BG,
                fontsize=9,
                fontweight="bold",
                va="center",
                fontfamily="monospace",
            )

            # ── Main plot axes ────────────────────────────────────────────────
            ax = fig.add_axes([0.07, 0.13, 0.92, 0.73])
            ax.set_facecolor(LCARS_SURFACE)

            # ── Bottom accent bar ─────────────────────────────────────────────
            ax_bot = fig.add_axes([0.07, 0.03, 0.92, 0.05])
            ax_bot.set_facecolor(LCARS_BG)
            ax_bot.axis("off")
            bot_segs = [
                (0.000, 0.040, LCARS_PURPLE),
                (0.042, 0.004, LCARS_BG),
                (0.046, 0.340, LCARS_ORANGE),
                (0.388, 0.004, LCARS_BG),
                (0.392, 0.180, LCARS_TEAL),
                (0.574, 0.004, LCARS_BG),
                (0.578, 0.422, LCARS_GOLD),
            ]
            for sx, sw, sc in bot_segs:
                ax_bot.add_patch(
                    plt.Rectangle(
                        (sx, 0),
                        sw,
                        1.0,
                        transform=ax_bot.transAxes,
                        color=sc,
                    )
                )

        # ── Solver loop ───────────────────────────────────────────────────────
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

        # ── Glow line helper ──────────────────────────────────────────────────
        def glow_plot(
            axes, x, y, color, lw=1.5, alpha_base=1.0, zorder=3, linestyle="-"
        ):
            for w, a in [(8, 0.04), (5, 0.09), (2.5, 0.2)]:
                axes.plot(
                    x,
                    y,
                    color=color,
                    lw=lw * w / 1.5,
                    alpha=alpha_base * a,
                    zorder=zorder - 1,
                    solid_capstyle="round",
                    linestyle="-",
                )
            axes.plot(
                x,
                y,
                color=color,
                lw=lw,
                alpha=alpha_base,
                zorder=zorder,
                solid_capstyle="round",
                linestyle=linestyle,
            )

        # ── Characteristics Family II ─────────────────────────────────────────
        for NII in range(1, int(n_max)):
            x_line, y_line = [0], [L]
            for NII_2 in range(1, NII + 1):
                x_NN, y_NN = grid.get_xy(NII - NII_2 + 1, NII_2)
                x_line.append(x_NN)
                y_line.append(y_NN)
                if NII == n_max - 1 and NII_2 != n_max - 1:
                    wall_x.append(x_NN)
                    wall_y.append(y_NN)
            if Graph2d:
                glow_plot(ax, x_line, y_line, LCARS_RED, lw=0.7, alpha_base=0.5)

        split_index = int(len(wall_x) * SP)
        wall_x = wall_x[:split_index]
        wall_y = wall_y[:split_index]

        # ── Characteristics Family I ──────────────────────────────────────────
        ks, ns = [], []
        for NI in range(1, n_max - 1):
            x_line, y_line = [], []
            for KI in range(1, int(k_max) - NI + 1):
                x_kn, y_kn = grid.get_xy(KI, NI)
                x_line.append(x_kn)
                y_line.append(y_kn)
            if 0 < GridField.get_x(grid, int(k_max) - NI, NI) < wall_x[-1]:
                if Graph2d:
                    glow_plot(ax, x_line, y_line, LCARS_PURPLE, lw=0.7, alpha_base=0.7)
            elif GridField.get_x(grid, int(k_max) - NI, NI) >= wall_x[-1]:
                if Graph2d:
                    glow_plot(ax, x_line, y_line, LCARS_TEAL, lw=0.7, alpha_base=0.5)
                ks.append(int(k_max) - NI + 1)
                ns.append(NI)

        y_calc = wall_y[-1]
        wall_x, wall_y = CC.CombustionChamber(wall_x, wall_y, R1, R2, L)
        wall_y_mirrored = -wall_y
        y_min = np.min(wall_y)

        A_calc = (y_calc / 1000) ** 2 * np.pi
        A_throat = (y_min / 1000) ** 2 * np.pi

        x_final_characteristic, _ = GridField.get_xy(grid, 1, n_max - 1)
        index = np.argmin(np.abs(wall_x - x_final_characteristic))
        y_final_characteristic = wall_y[index]

        A_exit = (y_final_characteristic / 1000) ** 2 * np.pi
        M_exit_characteristic = IT.AreaRatioInverse(A_calc / A_throat, g, "supersonic")

        P_exit = IT.Pressure(P_combustion, Param.g_exit, M_exit_characteristic)
        T_exit = IT.Temperature(T_combustion, Param.g_exit, M_exit_characteristic)
        A_sos = IT.LocalSoS(Param.g_exit, Rs, T_exit)

        Ve = A_sos * M_exit_characteristic * Efficiency
        Thrust = mdot * Ve + (P_exit - P_amb) * A_exit
        Exit_Angle = np.rad2deg(
            np.arctan2(wall_y[-1] - wall_y[-2], wall_x[-1] - wall_x[-2])
        )

        time1 = time.time()
        print(f"\nComputation Time: {time1 - time0:.2f} seconds\n")

        # ── 2D Plot ───────────────────────────────────────────────────────────
        if Graph2d:
            # Wall contour
            glow_plot(
                ax, wall_x, wall_y, LCARS_ORANGE, lw=2.2, alpha_base=1.0, zorder=5
            )
            glow_plot(
                ax,
                wall_x,
                wall_y_mirrored,
                LCARS_ORANGE,
                lw=2.2,
                alpha_base=1.0,
                zorder=5,
            )

            # Centreline dash
            x0 = np.array([np.min(wall_x), np.max(wall_x)])
            glow_plot(
                ax,
                x0,
                np.zeros(2),
                LCARS_GOLD,
                lw=0.9,
                alpha_base=0.7,
                zorder=4,
                linestyle="-.",
            )

            ax.set_xlabel("X  (mm)", fontsize=9)
            ax.set_ylabel("Y  (mm)", fontsize=9)
            ax.axis("equal")
            ax.grid(True)

            for spine in ax.spines.values():
                spine.set_edgecolor(LCARS_DIM)
                spine.set_linewidth(1.0)

            # Legend
            from matplotlib.lines import Line2D

            legend_elements = [
                Line2D([0], [0], color=LCARS_ORANGE, lw=2, label="NOZZLE WALL CONTOUR"),
                Line2D(
                    [0],
                    [0],
                    color=LCARS_RED,
                    lw=1,
                    linestyle="--",
                    alpha=0.8,
                    label="CHAR. FAMILY II",
                ),
                Line2D(
                    [0],
                    [0],
                    color=LCARS_PURPLE,
                    lw=1,
                    linestyle="--",
                    alpha=0.8,
                    label="CHAR. FAMILY I  (INTERNAL)",
                ),
                Line2D(
                    [0],
                    [0],
                    color=LCARS_TEAL,
                    lw=1,
                    linestyle="--",
                    alpha=0.8,
                    label="CHAR. FAMILY I  (EXIT)",
                ),
            ]
            leg = ax.legend(
                handles=legend_elements,
                loc="lower right",
                frameon=True,
                facecolor=LCARS_BG,
                edgecolor=LCARS_PURPLE,
                labelcolor=LCARS_GOLD,
                fontsize=8,
            )
            leg.get_frame().set_linewidth(1.5)

            # Data readout box

        # ── 3D Plot ───────────────────────────────────────────────────────────
        if Graph3d:
            fig = plt.figure(figsize=(16, 8))
            fig.set_facecolor(LCARS_BG)
            ax3 = fig.add_subplot(111, projection="3d")
            ax3.set_facecolor(LCARS_BG)

            theta = np.linspace(0, 2 * np.pi, 100)
            theta_grid, wall_x_grid = np.meshgrid(theta, wall_x)
            _, wall_y_grid = np.meshgrid(theta, wall_y)

            X = wall_x_grid
            Y = wall_y_grid * np.cos(theta_grid)
            Z = wall_y_grid * np.sin(theta_grid)

            # LCARS gradient: deep navy → purple → teal → gold
            lcars_cmap = LinearSegmentedColormap.from_list(
                "lcars",
                ["#0a0a14", "#9977cc", "#44aaaa", "#ffcc66"],
            )
            norm = plt.Normalize(wall_y.min() - 5, wall_y.max() + 5)
            colors = lcars_cmap(norm(wall_y_grid))
            colors[..., 3] = 0.92

            def update_plot(high_res=True):
                ax3.clear()
                ax3.set_facecolor(LCARS_BG)

                stride = (1, 1) if high_res else (8, 8)
                ax3.plot_surface(
                    X,
                    Y,
                    Z,
                    facecolors=colors,
                    linewidth=0,
                    antialiased=high_res,
                    rstride=stride[0],
                    cstride=stride[1],
                )

                # Profile lines in X-Z plane (both sides)
                ax3.plot(
                    wall_x,
                    np.zeros_like(wall_x),
                    wall_y,
                    color=LCARS_ORANGE,
                    lw=1.8,
                    alpha=0.95,
                    zorder=10,
                )
                ax3.plot(
                    wall_x,
                    np.zeros_like(wall_x),
                    wall_y,
                    color=LCARS_GOLD,
                    lw=5,
                    alpha=0.15,
                    zorder=9,
                )

                ax3.set_xlabel("X  (mm)", color=LCARS_GOLD, labelpad=8)
                ax3.set_ylabel("Y  (mm)", color=LCARS_GOLD, labelpad=8)
                ax3.set_zlabel("Z  (mm)", color=LCARS_GOLD, labelpad=8)
                ax3.tick_params(colors=LCARS_DIM, labelsize=7)

                for pane in (ax3.xaxis.pane, ax3.yaxis.pane, ax3.zaxis.pane):
                    pane.fill = True
                    pane.set_facecolor(LCARS_SURFACE)
                    pane.set_edgecolor(LCARS_DIM)

                ax3.set_title(
                    "NOZZLE GEOMETRY  //  3D SURFACE  //  MOC SOLVER",
                    color=LCARS_GOLD,
                    fontsize=11,
                    pad=14,
                    fontfamily="monospace",
                )

            x_range = wall_x.max() - wall_x.min()
            y_range = Y.max() - Y.min()
            z_range = Z.max() - Z.min()

            plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
            ax3.set_xlim(wall_x.min(), wall_x.max())
            ax3.set_ylim(Y.min(), Y.max())
            ax3.set_zlim(Z.min(), Z.max())
            ax3.set_box_aspect((x_range, y_range, z_range))

            def on_click(event):
                if event.inaxes == ax3:
                    update_plot(high_res=False)

            def on_release(event):
                update_plot(high_res=True)

            fig.canvas.mpl_connect("button_press_event", on_click)
            fig.canvas.mpl_connect("button_release_event", on_release)

            update_plot(high_res=True)

        print("\n" + "=" * 60)
        print("FINAL CALCULATION CHECK | SOLVER.PY")
        print("=" * 60)
        print(f"Input L (throat radius): {L:.4f} mm")
        print(f"Input mdot: {mdot:.6f} kg/s")
        print(f"Calculated y_min: {y_min:.4f} mm")
        print(f"A_throat: {A_throat * 1e6:.4f} mm²")
        print(f"A_exit: {A_exit * 1e6:.4f} mm²")
        print(f"Expansion ratio: {A_exit / A_throat:.4f}")
        print(f"M_exit: {M_exit_characteristic:.4f}")
        print(f"P_exit: {P_exit / 1e5:.4f} bar")
        print(f"T_exit: {T_exit:.2f} K")
        print(f"Ve: {Ve:.2f} m/s")
        print(f"Momentum thrust: {mdot * Ve:.2f} N")
        print(f"Pressure thrust: {(P_exit - P_amb) * A_exit:.2f} N")
        print(f"Thrust (before eff): {(mdot * Ve + (P_exit - P_amb) * A_exit):.2f} N")
        print(f"Efficiency: {Efficiency:.4f}")
        print(f"Thrust (after eff): {Thrust:.2f} N")
        print("=" * 60 + "\n")

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
            fig,
        )

    return solver()
