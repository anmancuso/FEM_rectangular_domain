# Finite Element Method for Homogeneous Poisson Equation

Python Scientific Developer Assessment - FEM 
This repository contains a Python implementation of the finite element method (FEM) for solving the homogeneous Poisson equation in a rectangular domain. The solver is implemented in solver.py, and its application is demonstrated in the Jupyter notebook FEM_Demo.ipynb.

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



