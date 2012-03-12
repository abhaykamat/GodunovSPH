#ifndef _GODUNOV1D_CPP
#define _GODUNOV1D_CPP

CGodunov1DSolverEulerBase::CGodunov1DSolverEulerBase() {
  time = 0.0;
}

CGodunov1DSolverEulerBase::~CGodunov1DSolverEulerBase() {}

// Compute the total energy, E, from primitive variables
// and ideal gas equation of state
void CGodunov1DSolverEulerBase::computeTotalEnergy(const double& d, const double& u, const double& p) {
  // See p88 of [1] (Toro, 2nd Edition).
  return d*(0.5*u*u + p/((gamma - 1)*d));
}

// Calculate the primitive pressure variable from the Total Energy, E
void CGodunov1DSolverEulerBase::computerPressureFromTotalEnergy(const double& d, const double& u, const double& E) {
  return (gamma - 1)*(E - 0.5*d*u*u);
}

// Specify precisely the left and right boundary conditions
bool CGodunov1DSolverEulerBase::computeDirichletBoundaryConditions(const double& dl, const double& dul, const double& El, const double& dr, const double& dur, const double& Er) {
  if (dl < 0.0 || El < 0.0 || dr < 0 || Er < 0) {   // Unphysical conditions
    return false;
  }
  
  // Specify left/right density
  d[0] = dl;
  d[cells-1] = dr;

  // Specify left/right fluid velocity
  du[0] = dul;
  du[cells-1] = dur;

  // Specify left/right fluid pressure
  E[0] = El;
  E[cells-1] = Er;
   
  return true;
}

void CGodunov1DSolverEulerBase::computePeriodicBoundaryConditions() {
  d[0] = d[cells-1];
  du[0] = du[cells-1];
  E[0] = E[cells-1];
}

void CGodunov1DSolverEulerBase::computeCFLCondition() {
  double smax = -1.0e6;
  
  // Find the maximum characteristic speed
  for (int i = 0; i < cells; i++) {
    if (fabs(du[i]/d[i]) > smax) {
      smax = fabs(du[i]/d[i])
    }
  }

  // Restrict the time step via CFL
  dt = cfl * dx / smax;

  // Check size of DT to avoid exceeding output time
  if ((time + dt) > timeout) {
    // Recompute DT
    dt = timeout - time
  }

}

// Update the global solution time step
void CGodunov1DSolverEulerBase::updateTimeStep() {
  time += dt;
}

void CGodunov1DSolverEulerBase::computeIntercellFluxes() {
  // Compute intercell flux for each respective quantity 
  // (density, velocity, pressure), at position i in 
  // vector flux_q where q in (d, du, E). Thus solution of
  // Riemann problem RP(i, i+1) is stored in flux_q(i).

  double ds, us, ps;
  
  for (int i = 0; i < cells; i++) {
    
    // Need to reformulate in primitive variable formulation 
    // to calculate Riemann problem solution
    double pl = computePressureFromTotalEnergy(d[i], u[i], p[i]);
    double pr = computePressureFromTotalEnergy(d[i+1], u[i+1], p[i+1]);

    // Solve the Riemann problem with the provided solver, making use of primitive variables
    CRiemannSolverEulerExact crs = CRiemannSolverEulerExact(d[i], du[i]/d[i], pl, d[i+1], du[i+1]/d[i+1], pr, gamma);
    crs.sampleWaveSolution(0.0, ds, us, ps);
    
    flux_d[i] = ds * us;
    flux_du[i] = ds * us * us + ps;
    flux_E[i] = us*(computeTotalEnergy(ds, us, ps) + ps);
    
  }

}

// Update the solution to a new time level using
// the explicit conservative formula 
void CGodunov1DSolverEulerBase::updateCells() {
  for (int i = 1; i < cells; i++) {
    d[i] += dtodx * (flux_d(i-1) - flux_d(i));
    du[i] += dtodx * (flux_du(i-1) - flux_du(i));
    E[i] += dtodx * (flux_E(i-1) - flux_E(i));
  }  
}

void CGodunov1DSolverEulerBase::outputCellData() {

}

void CGodunov1DSolverEulerBase::computeSolution() {
  for (int N = 1; N < N_max; n++) {

    // Set the Boundary Conditions
    if (bc = 0) {  // Periodic
      computePeriodicBoundaryConditions();
    } else {  // Dirichlet/fixed 
      computeDirichletBoundaryConditions();
    }
    
    // Algorithm update
    computeCFLCondition();       // Impose the Courant-Friedrichs-Lewy (CFL) condition
    updateTimeStep();            // Update the global time step
    computeIntercellFluxes();    // Compute intercell numerical fluxes
    updateCells();               // Update the solution according to conservative formula
    outputCellData();            // Determine whether to output the cell data and do so

  }  
}

#endif _GODUNOV1D_CPP

