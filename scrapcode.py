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
    #print(f"ve = {2 * (n_max-1) * dv}, Me = {IT.PMinv(2 * (n_max-1) * dv, g)}")

    
    print(f"[bold][purple]Test Area ratio: \t \t {A_test_ratio:.2f} \t \t | Test Mach: \t \t \t {Mach_test:.2f}[/bold][/purple]")
    print(f"[blue_violet]Test Thrust: \t \t \t {Thrust_test:.2f} N \t |")
    print(f"[bold][purple]> From Massflow: \t \t {mdot * Ve_test * Efficiency:.2f} N \t | > From Pressure: \t \t {Thrust_test - mdot * Ve_test * Efficiency:.2f} N\n")


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
    M_distribution = np.delete(M_distribution, -1)
    Exit_y = np.delete(Exit_y, -1)

    def MassWeightedThrust():
        coeffs = np.polyfit(Exit_y, M_distribution, 1)
        a, b = coeffs
        def M(y): return a * y + b

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