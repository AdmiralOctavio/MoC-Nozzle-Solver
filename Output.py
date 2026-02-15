from rich import print
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.rule import Rule
from rich.markup import escape


def m(s):
    return Text.from_markup(s)


def divider(color):
    return [Rule(style=color)]


import Parameters as Param
import matplotlib.pyplot as plt
import numpy as np
import stlgenerator
import os

folder_name = "Output_Files"
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

Graph = Param.Graph
Write = Param.Write
Stl = Param.Stl
Dxf = Param.Dxf

# Output file. Just outputs a table with the important design parameters, along with activating the write functions if required.

import solver


def outputTable():
    console = Console()
    (
        wall_x,
        Exit_Angle,
        y_min,
        T_combustion,
        g,
        P_combustion,
        M_optimal,
        AreaRatio,
        wall_y,
        M_exit,
        M_exit_characteristic,
        Thrust,
        A_exit,
        P_exit,
        Ve,
        ISP_cea,
        y_calc,
        L,
        wall_y_mirrored,
    ) = solver.solver(Graph)

    table = Table(
        show_header=False,
        box=None,
        pad_edge=True,
        expand=False,
        padding=(0, 2),
        collapse_padding=False,
    )
    table.add_column(max_width=50, justify="left")  # label 1
    table.add_column(max_width=50, justify="right")  # value 1
    table.add_column(min_width=5)  # SPACER COLUMN
    table.add_column(max_width=50, justify="left")  # label 2
    table.add_column(max_width=50, justify="right")  # value 2

    table.add_row(
        "[cyan3]Nozzle length:",
        f"[cyan3]{wall_x[-1]:.2f} mm",
        " ",
        "[cyan3]Total length:",
        f"[cyan3]{wall_x[-1] - wall_x[0]:.2f} mm",
    )

    exit_style = "red" if Exit_Angle > 6 else "dark_turquoise"
    table.add_row(
        m("[dark_turquoise]Exit Angle:[/dark_turquoise]"),
        m(f"[{exit_style}]{Exit_Angle:.2f} deg"),
        " ",
        m("[dark_turquoise]Exit radius:"),
        f"[dark_turquoise]{wall_y[-1]:.2f} mm",
    )

    table.add_row(
        "[cyan3]True throat radius:",
        f"[cyan3]{y_min:.2f} mm",
        " ",
        "[cyan3]True throat diameter:",
        f"[cyan3]{2*y_min:.2f} mm",
    )

    table.add_row(
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
    )

    table.add_row(
        "[orange_red1]Combustion temperature:",
        f"[orange_red1]{T_combustion:.1f} K",
        " ",
        "[orange_red1]Gamma:",
        f"[orange_red1]{g:.3f}",
    )

    table.add_row(
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
    )

    table.add_row(
        "[light_sky_blue3]Optimal pressure ratio:",
        f"[light_sky_blue3]{P_combustion/101325:.2f}",
        "",
        "",
    )
    table.add_row(
        "[sky_blue2]Optimal exit Mach:", f"[sky_blue2]{M_optimal:.2f}", "", ""
    )
    table.add_row(
        "[light_sky_blue3]Optimal expansion ratio:",
        f"[light_sky_blue3]{AreaRatio:.2f}",
        "",
        "",
    )

    table.add_row(
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
        Rule(style="cyan"),
    )

    table.add_row(
        "[green_yellow]Theoretical expansion ratio:",
        f"[green_yellow]{(wall_y[-1]**2 / L**2):.2f}",
        " ",
        "[green_yellow]True expansion ratio:",
        f"[green_yellow]{(y_calc**2 / y_min**2):.2f}",
    )

    table.add_row(
        "[light_green]Design exit Mach:",
        f"[light_green]{M_exit}",
        " ",
        "[light_green]Predicted exit Mach:",
        f"[light_green]{M_exit_characteristic:.2f}",
    )

    table.add_row(
        "[green_yellow]Predicted thrust:", f"[green_yellow]{Thrust:.0f} N", "", ""
    )

    table.add_row(
        "[light_green]> Thrust (massflow):",
        f"[light_green]{Thrust - (P_exit-101325)*A_exit:.2f} N",
        " ",
        "[light_green]> Thrust (pressure):",
        f"[light_green]{(P_exit-101325)*A_exit:.2f} N",
    )

    if P_exit < 0.25 * 101325:
        warn = Text("FLOW SEPARATION WILL OCCUR", style="bold red")
    elif P_exit < 0.4 * 101325:
        warn = Text("Flow separation may occur", style="bold orange1")
    elif P_exit < 0.5 * 101325:
        warn = Text("Nearing exit instability region", style="yellow")
    else:
        warn = Text("Flow stable", style="green")

    table.add_row(
        "[green_yellow]Exit pressure:", f"[green_yellow]{P_exit:.0f} Pa", "", ""
    )
    table.add_row(warn, "", "", "")

    table.add_row(
        "[light_green]Specific impulse:",
        f"[light_green]{Ve/9.80665:.2f} s",
        " ",
        "[light_green]CEA ISp:",
        f"[light_green]{ISP_cea:.2f} s",
    )

    panel = Panel(
        table,
        title="[bold red]Output Nozzle Design Specifications[/bold red]",
        border_style="cyan",
        padding=(1, 2),
        expand=True,
    )

    console.print(panel)

    if Write == True:
        combined = np.stack((wall_x / 1000, wall_y / 1000), axis=1)
        filename = f"Nozzle_Contour_M={M_exit:.2f}.csv"
        filepath = os.path.join(folder_name, filename)
        np.savetxt(
            filepath,
            combined,
            delimiter=",",
            header="x_position (mm), y_radius (mm)",
            comments="",
        )

    if Stl == True:
        stlgenerator.create_stl(wall_x, wall_y, M_exit)

    if Dxf == True:
        stlgenerator.create_dxf(wall_x / 1000, wall_y_mirrored / 1000, M_exit)
