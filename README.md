# Finite Element Method for Homogeneous Poisson Equation

### Python Scientific Developer Assessment - FEM 

This repository contains a Python implementation of the finite element method (FEM) for solving the homogeneous Poisson equation in a rectangular domain. The solver is implemented in `solver.py`, and its application is demonstrated in the Jupyter notebook `FEM_Demo.ipynb`.

Description
The solver addresses the following problem:

$$
-Δu(x) = f(x) \quad \forall x \in \Omega := [a, b] \times [c, d]
$$

with the boundary condition:

$$
u(x) = 0 \quad \forall x \in \partial\Omega
$$


The implementation focuses on:

- Using essential libraries like scipy and numpy for computations.
- Analyzing the method's convergence with a test case.
- Visualizing results using matplotlib or plotly.
- Refining meshes based on convergence theory to enhance solution accuracy.

## Files in the Repository
- `solver.py`: Contains all the key functions for the FEM implementation.
- `FEM_Demo.ipynb`: A Jupyter notebook demonstrating the usage of functions from `solver.py`.
- 
## Usage
To use the solver, run the Jupyter notebook `FEM_Demo.ipynb` . This notebook will guide you through the process of solving the Poisson equation using the functions defined in `solver.py`.

## Installation
To run the code, ensure you have Python installed along with the following libraries:

```
numpy
scipy
matplotlib
```

You can install these dependencies using pip:
```
 pip install numpy scipy matplotlib
```

## Key Functions
### solver.py


#### `create_delaunay_mesh(a, b, c, d, nx, ny)`
Generates a Delaunay triangulation mesh over a rectangular grid.
- **Parameters**: `a, b` (x-axis range), `c, d` (y-axis range), `nx, ny` (divisions along the x and y axes).
- **Returns**: `points` (2D array of points in the grid), `delaunay_tri.simplices` (the simplices, i.e., triangles, of the Delaunay triangulation).

#### `plot_delaunay_mesh(points, triangles)`
Plots the Delaunay triangulation mesh.
- **Parameters**: `points` (2D array of points in the grid), `triangles` (the simplices of the Delaunay triangulation).

#### `find_boundary_nodes_delaunay(points, a, b, c, d)`
Identifies boundary nodes in a set of points within a rectangular region.
- **Parameters**: `points` (2D array of points), `a, b` (x-coordinate boundaries), `c, d` (y-coordinate boundaries).
- **Returns**: A sorted list of indices of points that lie on the boundary of the rectangle.

#### `compute_local_stiffness_matrix(coords)`
Calculates the local stiffness matrix for a triangular element. The local stiffness matrix $`K`$ for a triangle is calculated using the coordinates of its vertices. The general formula is:

$$
K_{ij} = \int_{\Omega} \nabla \phi_i \cdot \nabla \phi_j , d\Omega
$$

where $`phi_i`$ and $`\phi_j`$ are the linear shape functions associated with the triangle's nodes. The area $`\Omega`$ of the triangle is given by:

$$
\text{Area} = \frac{1}{2} |x_1(y_2 - y_3) + x_2(y_3 - y_1) + x_3(y_1 - y_2)|
$$

- **Parameters**: `points` (3x2 array with the (x, y) coordinates of the triangle's vertices).
- **Returns**: 3x3 local stiffness matrix.

#### `transform_stiffness_matrix(points, K_ref)`
Transforms a reference stiffness matrix into the stiffness matrix for a generic triangle. This transformation is necessary to adapt the stiffness matrix calculated for a reference triangle to any triangle in the domain. The transformation involves an affine mapping between the reference triangle and the generic triangle, which is described by a linear transformation and a translation.

The affine transformation is defined as follows:
- Let the vertices of the reference triangle be `A_ref`, `B_ref`, and `C_ref`, and the vertices of the generic triangle be `A`, `B`, and `C`.
- The transformation maps the vertices of the reference triangle to those of the generic triangle.
- This mapping can be represented by a matrix `T` such that `T * [A_ref; B_ref; C_ref] = [A; B; C]`.

The Jacobian $`J`$ of the transformation, which is the matrix of all first-order partial derivatives, is used to transform the area element and the gradient of the shape functions. The determinant of the Jacobian $`\det(J)`$ represents the area scaling factor between the reference and the generic triangle.

The transformed stiffness matrix `K` for the generic triangle is then calculated using the formula:

$$
K = det(J) * T^{-1} ^T * K_{ref} * T^{-1}
$$

where $`K_{ref}`$ is the stiffness matrix for the reference triangle.

- **Parameters**: `points` (coordinates of the generic triangle), `K_ref` (reference stiffness matrix).
- **Returns**: Transformed stiffness matrix for the generic triangle.

#### `assemble_global_stiffness_matrix(triangles, nodes, K_ref)`
Assembles the global stiffness matrix from the local stiffness matrices of triangles. The global stiffness matrix is a sparse matrix representing the entire FEM system. Each local stiffness matrix contributes to the global matrix at positions corresponding to the triangle's nodes.

- **Parameters**: `triangles` (list of triangles), `nodes` (node coordinates), `K_ref` (reference stiffness matrix).
- **Returns**: `lil_matrix`: Global stiffness matrix (expected to be a sparse matrix).
### FEM_Demo.ipynb








## Contributing
Contributions to enhance the solver's efficiency, accuracy, or functionality are welcome. Please feel free to fork the repository and submit pull requests.
