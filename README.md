# GodunovSPH

GodunovSPH is a Computational Fluid Dynamics library for modelling compressible flows in physics and engineering. In particular, it enables numerical experimentation in aeronautics, astronautics and astrophysics. It is based on a Generalsed Riemann Solver/Godunov Finite Volume Algorithm for shocks as well as Smooth Particle Hydronamics for astrophysical flow modelling.

## Languages

* The core solver routines will be written in C++
* Depending upon complexity, pre- and post-analysis data manipulation will be written in Python

## Roadmap

* Create CRiemannSolverEulerBase abstract base class to be inherited by all subsequent Riemann solver objects, based on Euler equations of gas dynamics
* Create an exact Riemann solver (CRiemannSolverEulerExact) for the one-dimensional Euler equations based on Chapter 4 of [1], providing a testing capability
* Create a Harten, Lax and van Leer Riemann solver (HLL) as well as Harten, Lax and van Leer with contact capability (via Toro); HLLC. Finally, create the HLLE solver (Einfeldt modification)
* Create a Riemann solver based on Roe's approximation, which will be useful when creating semi-homogeneous flow fields
* Create a first order one-dimensional Godunov solver for the Euler equations that utilises a Riemann solver to determine the inter-cell fluxes, based on Chapter 3, 4, 5 and 6 of [1]
* Create a second order Riemann solver based on the Generalised Riemann Problem (GRP) of Ben-Artzi and Falcovitz
* Create a second order (space/time) Godunov solver for the 1D Euler equations using the GRP described above
* Extend the 1D GRP-based Godunov solver to multiple dimensions (2D to begin with) on Cartesian structured geometries
* Extend the GRP-based Godunov solver to unstructured meshes, including a variety of test problems
* Create a basic smoothed particle hydrodynamics (SPH) code, source TBD
* Add the GRP Riemann solver to the SPH code to handle flow discontinuities
* Apply the software to multiple test objects in aeronautics, astronautics and astrophysics

## References

[1] Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction, 2nd Edition, E.F. Toro