# Diffusion_Equation
Finite Difference Methods for Solving the Diffusion Equation

# Numerical Analysis of Diffusion Equation with Finite Difference Methods

This project presents a numerical solution to the diffusion equation using finite difference methods in C++. It explores the behavior of the diffusion equation for various values of the time-stepping parameter \( \theta \) and stability parameter \( r \). The project also includes Python scripts to visualize the results.

## Overview

The diffusion equation is a key partial differential equation in fluid dynamics, thermodynamics, and other fields. Solving it numerically allows for an understanding of diffusion behavior over time and space. This project implements three common finite difference methods:
1. **Explicit Scheme ( \( \theta = 0 \))**
2. **Crank-Nicholson Scheme ( \( \theta = 0.5 \))**
3. **Implicit Scheme ( \( \theta = 1 \))**

For each method, the parameter \( r \), which represents the ratio of time step size to the spatial step size squared, is varied.

## Features
- Numerical solution for the diffusion equation in C++.
- Support for multiple values of \( \theta \) and \( r \).
- Error history tracking and solution history visualization.
- Python scripts for plotting error and solution histories.

## Files and Structure

- **Numerical_Analyis.cpp**: The main C++ file that implements the finite difference methods for solving the diffusion equation.
- **solution_r_X_theta_Y.txt**: Solution files that store \( u(x,t) \) values for different combinations of \( r \) and \( \theta \).
- **error_r_X_theta_Y.txt**: Error history files that store the error (L-inf norm) over time for each combination of \( r \) and \( \theta \).
- **plot_error_solution.py**: Python script that reads the solution and error files and generates plots for visualization.
- **README.md**: This file, provides an overview of the project.

## Usage

### Running the C++ Code
1. Compile the C++ code using a C++ compiler like `g++`:
   ```bash
   g++ Numerical_Analyis.cpp -o numerical_analysis

