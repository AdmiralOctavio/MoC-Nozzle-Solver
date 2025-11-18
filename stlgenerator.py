import numpy as np
from stl import mesh
import ezdxf

def create_dxf(x_contour, y_contour, Mach):
    filename = f"Nozzle_Contour_M={Mach:.2f}.dxf"
    doc = ezdxf.new(dxfversion='R2018')
    msp = doc.modelspace()

    points = []
    for x, y in zip(x_contour, y_contour):
        points.append((x, y, 0))

    msp.add_lwpolyline(points, 
                       dxfattribs={
                           'layer': 'NOZZLE_CONTOUR',
                           'color': 256 # ByLayer color
                       })
    
    # 4. Save the file
    doc.saveas(filename)


def create_stl(x_contour, y_contour, Mach):
    filename = f"Nozzle_Contour_M={Mach:.2f}.stl"
    num_points = len(x_contour)
    num_facets = 60
    theta = np.linspace(0, 2 * np.pi, num_facets, endpoint=False) 

    num_triangles = 2 * (num_points - 1) * num_facets
    
    nozzle_mesh = np.zeros(num_triangles, dtype=mesh.Mesh.dtype)

    vertices = np.zeros((num_points, num_facets, 3))
    
    for i in range(num_points):
        x_val = x_contour[i]
        r_val = y_contour[i]
        
        vertices[i, :, 0] = x_val
        vertices[i, :, 1] = r_val * np.cos(theta)
        vertices[i, :, 2] = r_val * np.sin(theta)

    triangle_index = 0
    
    for i in range(num_points - 1):
        for j in range(num_facets):
            v1 = vertices[i, j]
            v2 = vertices[i, (j + 1) % num_facets] 
            
            v3 = vertices[i + 1, (j + 1) % num_facets]
            v4 = vertices[i + 1, j]

            nozzle_mesh['vectors'][triangle_index] = [v1, v2, v3]
            triangle_index += 1
            
            nozzle_mesh['vectors'][triangle_index] = [v1, v3, v4]
            triangle_index += 1

    mesh_obj = mesh.Mesh(nozzle_mesh)
        
    mesh_obj.save(filename)