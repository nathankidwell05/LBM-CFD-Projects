# Lid-Driven Cavity Flow

## Description

This project implements a D2Q9 lattice Boltzmann solver for the lid-driven cavity flow benchmark at Reynolds number 100.

Features:
- D2Q9 lattice
- MRT collision operator with selectable relaxation cases
- Moving-lid boundary condition
- No-slip stationary walls
- Streamline, velocity, density, and Ghia benchmark comparison plots

## Governing Equations

Macroscopic density:

rho = sum(f_i)

Macroscopic velocity:

u = sum(f_i e_i) / rho

MRT collision:

f* = M^-1 [m - S(m - m_eq)]

## Files

- CavityFlow.m
- eqm_d2q9.m
- moment_rho_u_d2q9.m
- plot_streamlines.m

## Skills Demonstrated

- MATLAB
- Computational Fluid Dynamics
- Lattice Boltzmann Method
- MRT Collision Modeling
- Boundary Conditions
- Benchmark Validation
