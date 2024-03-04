import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse import lil_matrix
import matplotlib.pyplot as plt
import matplotlib.patches as patches





def create_delaunay_mesh(a, b, c, d, nx, ny):
    """
    Creates a Delaunay triangulation mesh over a rectangular grid.

    Parameters:
    a, b - Define the range of x-axis values.
    c, d - Define the range of y-axis values.
    nx, ny - Define the number of divisions along the x and y axes.

    Returns:
    points - A 2D array of points in the mesh grid.
    delaunay_tri.simplices - The simplices (triangles) of the Delaunay triangulation.
    """
    
    x = np.linspace(a, b, nx)
    y = np.linspace(c, d, ny)
    points = np.array(np.meshgrid(x, y)).T.reshape(-1, 2)# 2D array of (x, y) from the meshgrid of x and y values.
    delaunay_tri = Delaunay(points)  
    #Delaunay da scipy : usato in ambiti di grafica computerizzata, simulazioni numeriche e GIS.

    return points, delaunay_tri.simplices


def plot_delaunay_mesh(points, triangles):
    """
    Plots the Delaunay triangulation mesh.

    Parameters:
    points - A 2D array of points in the mesh grid.
    triangles - The simplices (triangles) of the Delaunay triangulation.
    
    This function creates a plot of the Delaunay triangulation mesh using matplotlib.
    """
    x_min, y_min = points.min(axis=0)
    x_max, y_max = points.max(axis=0)
    width = x_max - x_min
    height = y_max - y_min
    fig, ax = plt.subplots()
    
    ax.set_xlim(x_min - width * 0.1, x_max + width * 0.1)
    ax.set_ylim(y_min - height * 0.1, y_max + height * 0.1)
    plt.triplot(points[:, 0], points[:, 1], triangles, '.-', label="Nodes")
    rect = patches.Rectangle((x_min, y_min), width, height, linewidth=2, edgecolor='r', facecolor='none', label="Border")
    ax.add_patch(rect)    
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1))
    plt.show()


def find_boundary_nodes_delaunay(points, a, b, c, d):
    """
    Identifies boundary nodes in a set of points within a rectangular region. This will be used to apply the boundary conditions.

    Parameters:
    points - A 2D array of points.
    a, b - The x-coordinate boundaries of the rectangle.
    c, d - The y-coordinate boundaries of the rectangle.

    Returns:
    A sorted list of indices of points that lie on the boundary of the rectangle.
    """
    
    boundary_nodes = set()  

    # Iterate over each point and its index
    for i, point in enumerate(points):
        x, y = point  

        
        if x == a or x == b or y == c or y == d:
            boundary_nodes.add(i) 

    return sorted(boundary_nodes)  



def compute_local_stiffness_matrix(coords):
    """
    Compute the local stiffness matrix for a triangular element; this is used to compute the local matrix for a simple basic triangle.
    
    
    Parameters:
    points (np.array): 3x2 array containing the (x, y) coordinates of the triangle's vertices.

    Returns:
    np.array: 3x3 local stiffness matrix.
    """
    
    x1, y1 = points[0]
    x2, y2 = points[1]
    x3, y3 = points[2]
    
    
    area = 0.5 * abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))

    # Gradienti funzioni di forma
    b = np.array([y2 - y3, y3 - y1, y1 - y2]) / (2 * area)
    c = np.array([x3 - x2, x1 - x3, x2 - x1]) / (2 * area)

    # Matrice di rigidità
    K = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            K[i, j] = area * (b[i] * b[j] + c[i] * c[j])

    return K



def transform_stiffness_matrix(points, K_ref):
    """
    Transform a reference stiffness matrix into the stiffness matrix for a generic triangle.

    Parameters:
    points (np.array): 3x2 array containing the (x, y) coordinates of the generic triangle's vertices.
    K_ref (np.array): Reference stiffness matrix.

    Returns:
    np.array: Transformed stiffness matrix for the generic triangle.
    """
    # Coordinates of the reference triangle
    ref_points = np.array([[0, 0], [1, 0], [0, 1]])
    
    # Area of the reference triangle
    area_ref = 0.5
    
    # Area of the generic triangle
    x1, y1 = points[0]
    x2, y2 = points[1]
    x3, y3 = points[2]
    area = 0.5 * abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
    
    # Compute the transformation matrix B and vector b
    A = np.column_stack([points[:, 0], points[:, 1], np.ones(3)])
    B = np.linalg.inv(A) @ ref_points

    # Jacobian of the transformation
    J = B[:2, :2]
    det_J = np.linalg.det(J)

    # Scale factor for the area transformation
    scale = area / area_ref

    # Transform the stiffness matrix
    K_transformed = scale * (np.linalg.inv(J).T @ K_ref[:2, :2] @ np.linalg.inv(J))

    # Expand the transformed stiffness matrix to 3x3
    K_transformed_full = np.zeros_like(K_ref)
    K_transformed_full[:2, :2] = K_transformed
    K_transformed_full[2, 2] = scale * K_ref[2, 2]  # Handle the (3,3) entry separately

    return K_transformed_full





def assemble_global_stiffness_matrix(triangles, nodes, K_ref):
    """
    Assemble the global stiffness matrix from local stiffness matrices of triangles.

    Parameters:
    triangles (list of lists): List of triangles, each represented by a list of node indices.
    nodes (np.array): Array of node coordinates.
    K_ref (np.array): Reference stiffness matrix.

    Returns:
    lil_matrix: Global stiffness matrix. Since from the theory we expect it to be a sparse matrix.
    """
    n_nodes = len(nodes)
    K_global = lil_matrix((n_nodes, n_nodes))

    for tri in triangles:
        points = np.array([nodes[tri[0]], nodes[tri[1]], nodes[tri[2]]])

        # Calculate the local stiffness matrix for the current triangle
        K_local = transform_stiffness_matrix(points, K_ref)

        # Add the local stiffness matrix to the global matrix
        for i in range(3):
            for j in range(3):
                K_global[tri[i], tri[j]] += K_local[i, j]

    return K_global




