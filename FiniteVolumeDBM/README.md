# Finite-Volume Discrete Boltzmann Method

## Description

This project implements a D2Q9 finite-volume discrete Boltzmann method (FVDBM) for the lid-driven cavity flow benchmark at Reynolds number 100. The MATLAB solver uses cell-centered control volumes, first-order upwind face fluxes, BGK collision, and non-equilibrium extrapolation at the cavity walls.

Features:
- D2Q9 discrete velocity model
- Cell-centered finite-volume discretization
- First-order upwind face fluxes
- BGK collision model
- Moving-lid and no-slip wall boundary conditions
- Low-Mach and explicit time-step stability checks
- Ghia centerline benchmark comparison
- Automatic velocity, density, and convergence plots

## Governing Equation

The discrete distribution functions satisfy:

df_i/dt + xi_i . grad(f_i) = -(f_i - f_i^eq) / tau

Macroscopic density and velocity are recovered from:

rho = sum(f_i)

u = sum(f_i xi_i) / rho

## Numerical Setup

The default case uses:

- Reynolds number: `Re = 100`
- Grid: `81 x 81` control volumes
- Lid velocity: `U_lid = 0.1`
- Maximum iterations: `100000`
- Convergence tolerance: `1e-8`

The saved case converged after 61,000 iterations with a final velocity residual of approximately `7.2e-9`. Mean absolute centerline errors relative to Ghia et al. are `0.0217` for `u/U_lid` and `0.0231` for `v/U_lid`.

## Results

### Velocity Magnitude

![Velocity Magnitude](velocity_magnitude.png)

### Velocity Field

![Velocity Field](velocity_vectors.png)

### Density Field

![Density Field](density_field.png)

### Centerline Benchmark Comparison

![Centerline Benchmark Comparison](centerline_comparison.png)

### Convergence History

![Convergence History](convergence_history.png)

## Files

- `FVBDM_v2.m`
- `FVBDM_results.mat`
- `velocity_magnitude.png`
- `velocity_vectors.png`
- `density_field.png`
- `centerline_comparison.png`
- `convergence_history.png`

## Reference

U. Ghia, K. N. Ghia, and C. T. Shin, "High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method," Journal of Computational Physics, vol. 48, no. 3, pp. 387-411, 1982.

## Skills Demonstrated

- MATLAB
- Computational Fluid Dynamics
- Finite-Volume Method
- Discrete Boltzmann Method
- BGK Collision Modeling
- Boundary Conditions
- Benchmark Validation
