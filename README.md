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




