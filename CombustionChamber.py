import numpy as np
import Parameters as P

Length = P.L_combustion
Diameter = P.D_combustion

def CombustionChamber(wall_x, wall_y, x_arc, y_arc, R, n_fillet=20):

    def chamber_slope(y):
        return -1 * (y - y_arc[0]) + x_arc[0]

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
    P1_shift = P1 + R * n1
    
    c1p = -(a1 * P1_shift[0] + b1 * P1_shift[1])

    P2 = np.array([x_corner, y_corner])
    P2_shift = P2 + R * n2
    c2p = -(a2 * P2_shift[0] + b2 * P2_shift[1])

    A = np.array([[a1, b1],
                  [a2, b2]])
    C = -np.array([c1p, c2p])
    
    x_center, y_center = np.linalg.solve(A, C)
    center = np.array([x_center, y_center])

    
    def project_point_to_line(point, a, b, c):
        
        denom = a*a + b*b
        x0, y0 = point
        t = (a*x0 + b*y0 + c) / denom
        
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

    x_fillet = x_center + R * np.cos(theta)
    y_fillet = y_center + R * np.sin(theta)

    y_contract = np.linspace(tangent_s[1], y_arc[0], 2)
    x_contract = chamber_slope(y_contract)

    chamber_x = np.linspace(x_fillet[0] - 100.0, x_fillet[0], 2)
    chamber_y = np.ones(2) * y_fillet[0]

    
    wall_x = np.concatenate([chamber_x, x_fillet, x_contract, wall_x])
    wall_y = np.concatenate([chamber_y, y_fillet, y_contract, wall_y])

    return wall_x, wall_y
