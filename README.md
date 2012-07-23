# GodunovSPH

GodunovSPH is a Computational Fluid Dynamics library for modelling compressible flows in physics and engineering. In particular, it enables numerical experimentation in aeronautics, astronautics and astrophysics. 

It is based on a Generalsed Riemann Problem quasi-conservative scheme as an extension of Godunov's Solver (Finite Volume Method) for discontinuity capturing (shocks and contacts).

It also enables astrophysical modelling via a Smoothed Particle Hydrodynamics code, which makes use of the GRP solver to capture hyperbolic effects.

## Code

* The core solver routines are written in object-oriented C++

## Roadmap

* Create RiemannSolverEulerBase abstract base class to be inherited by all subsequent Riemann solver objects, based on Euler equations of gas dynamics
* Create an exact Riemann solver (RiemannSolverEulerExact) for the one-dimensional Euler equations based on Chapter 4 of [1], providing a testing capability
* Create a second order (space/time) Godunov-type conservative scheme solver for the 1D Euler equations using the Generalised Riemann Problem method of Ben-Artzi and Falcovitz [2]
* Extend the 1D GRP-based Godunov solver to multiple dimensions (2D to begin with) on Cartesian structured geometries, based on the operator splitting techniques of Strang
* Extend the GRP-based Godunov solver to unstructured meshes, including a variety of test problems
* Create a smoothed particle hydrodynamics (SPH) code
* Add the GRP Riemann solver to the SPH code to handle flow discontinuities
* Apply the software to multiple test objects in aeronautics, astronautics and astrophysics

## Extensions

* Create a Harten, Lax and van Leer Riemann solver (HLL) as well as Harten, Lax and van Leer with contact capability (via Toro); HLLC. Finally, create the HLLE solver (Einfeldt modification)
* Create a Riemann solver based on Roe's approximation, which will be useful when creating semi-homogeneous flow fields

## References

[1] Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction, 2nd Edition, E.F. Toro
[2] Generalized Riemann Problems in Computational Fluid Dynamics, M. Ben-Artzi and J. Falcovitz