import numpy as np
import Parameters as P

Length = P.L_combustion
L_throat = P.L
alpha = np.deg2rad(P.Chamber_Slope)


def CombustionChamber(wall_x, wall_y, R1, R2):

    Rt = L_throat 
    Contraction_Ratio = P.Contraction_ratio
    R_inlet = Rt * np.sqrt(Contraction_Ratio) 

    dy_straight = (R_inlet - Rt) - (R1 + R2) * (1 - np.cos(alpha))
    L_straight = dy_straight / np.sin(alpha)

    L_total = (R1 + R2) * np.sin(alpha) + L_straight * np.cos(alpha)

    theta_inlet = np.linspace(0, alpha, 50)
    x_arc_inlet = -L_total + R2 * np.sin(theta_inlet)
    y_arc_inlet = R_inlet - R2 * (1 - np.cos(theta_inlet))

    x_line_cont = -Length
    y_line_cont = R_inlet

    x_t1, y_t1 = x_arc_inlet[-1], y_arc_inlet[-1]
    x_t2 = x_t1 + L_straight * np.cos(alpha)
    y_t2 = y_t1 - L_straight * np.sin(alpha)
    x_line = np.array([x_t1, x_t2])
    y_line = np.array([y_t1, y_t2])

    theta_throat = np.linspace(alpha, 0, 50)
    x_arc_throat = -R1 * np.sin(theta_throat)
    y_arc_throat = Rt + R1 * (1 - np.cos(theta_throat))

    wall_x = np.append(x_arc_throat, wall_x)
    wall_x = np.append(x_line, wall_x)
    wall_x = np.append(x_arc_inlet, wall_x)
    wall_x = np.append(x_line_cont, wall_x)

    wall_y = np.append(y_arc_throat, wall_y)
    wall_y = np.append(y_line, wall_y)
    wall_y = np.append(y_arc_inlet, wall_y)
    wall_y = np.append(y_line_cont, wall_y)

    return wall_x, wall_y
