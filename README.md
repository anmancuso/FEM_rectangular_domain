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



## Key Functions
### solver.py

### FEM_Demo.ipynb


## Installation
To run the code, ensure you have Python installed along with the following libraries:

```
numpy
scipy
matplotlib (or plotly for advanced visualizations)
```

You can install these dependencies using pip:
```
 pip install numpy scipy matplotlib
```


## Usage
To use the solver, run the Jupyter notebook `FEM_Demo.ipynb` . This notebook will guide you through the process of solving the Poisson equation using the functions defined in `solver.py`.


## Contributing
Contributions to enhance the solver's efficiency, accuracy, or functionality are welcome. Please feel free to fork the repository and submit pull requests.
