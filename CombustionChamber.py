import numpy as np
import Parameters as P
import matplotlib.pyplot as plt

Length = P.L_combustion
Diameter = P.D_combustion
L_throat = P.L
alpha = np.deg2rad(P.Chamber_Slope)


def CombustionChamber(wall_x, wall_y, R1, R2):
    n_fillet = 20
    wall_x = np.array(wall_x)
    wall_y = np.array(wall_y)
    wall_x = np.append(0, wall_x)
    wall_y = np.append(L_throat, wall_y)

    # First arc immediately before the nozzle
    x_arc1 = np.linspace(-R1 * np.sin(alpha), wall_x[0], 20)
    y_arc1 = -np.sqrt(R1**2 - x_arc1**2) + L_throat + R1

    wall_x = np.append(x_arc1, wall_x)
    wall_y = np.append(y_arc1, wall_y)

    # Second, straight arc between nozzle and first fillet
    x_arc2 = np.linspace(-Diameter / 2 * np.tan(alpha), x_arc1[0], 5)
    y_arc2 = np.linspace(Diameter / 2, y_arc1[0], 100)

    # Third straight arc, just combustion chamber contour
    x_arc3 = np.linspace(-Length - x_arc2[0], x_arc2[0], 20)
    y_arc3 = np.linspace(Diameter / 2, y_arc2[0], 100)

    def chamber_slope(y):
        return -np.tan(alpha) * (y - y_arc1[0]) + x_arc1[0]

    Upper_bound = Diameter / 2

    x_corner = chamber_slope(Upper_bound)
    y_corner = Upper_bound

    a1, b1, c1 = 0.0, 1.0, -Upper_bound

    a2, b2, c2 = 1.0, 1.0, -(x_corner + y_corner)

    n1 = np.array([0.0, -1.0])

    n2_candidate = np.array([a2, b2], dtype=float)
    n2_unit = n2_candidate / np.linalg.norm(n2_candidate)

    if np.dot(n2_unit, n1) < 0:
        n2_unit = -n2_unit
    n2 = n2_unit

    P1 = np.array([0.0, Upper_bound])
    P1_shift = P1 + R2 * n1

    c1p = -(a1 * P1_shift[0] + b1 * P1_shift[1])

    P2 = np.array([x_corner, y_corner])
    P2_shift = P2 + R2 * n2
    c2p = -(a2 * P2_shift[0] + b2 * P2_shift[1])

    A = np.array([[a1, b1], [a2, b2]])
    C = -np.array([c1p, c2p])

    x_center, y_center = np.linalg.solve(A, C)
    center = np.array([x_center, y_center])

    def project_point_to_line(point, a, b, c):

        denom = a * a + b * b
        x0, y0 = point
        t = (a * x0 + b * y0 + c) / denom

        xf = x0 - a * t
        yf = y0 - b * t
        return np.array([xf, yf])

    tangent_h = project_point_to_line(center, a1, b1, c1)
    tangent_s = project_point_to_line(center, a2, b2, c2)

    d1 = np.linalg.norm(center - tangent_h)
    d2 = np.linalg.norm(center - tangent_s)

    ang_h = np.arctan2(tangent_h[1] - y_center, tangent_h[0] - x_center)
    ang_s = np.arctan2(tangent_s[1] - y_center, tangent_s[0] - x_center)

    def angle_diff(a, b):
        diff = b - a
        diff = (diff + np.pi) % (2 * np.pi) - np.pi
        return diff

    diff = angle_diff(ang_h, ang_s)
    theta = np.linspace(ang_h, ang_h + diff, n_fillet)

    x_fillet = x_center + R2 * np.cos(theta)
    y_fillet = y_center + R2 * np.sin(theta)

    y_contract = np.linspace(tangent_s[1], wall_y[0], 2)
    x_contract = chamber_slope(y_contract)

    chamber_x = np.linspace(x_fillet[0] - Length, x_fillet[0], 2)
    chamber_y = np.ones(2) * y_fillet[0]

    wall_x = np.concatenate([chamber_x, x_fillet, x_contract, wall_x])
    wall_y = np.concatenate([chamber_y, y_fillet, y_contract, wall_y])

    return wall_x, wall_y
