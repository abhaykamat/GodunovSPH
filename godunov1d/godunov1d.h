#ifndef _GODUNOV1D_H
#define _GODUNOV1D_H

#include "riemann.h"

using namespace std;

class CGodunov1DSolverEulerBase {
 protected:
  double gamma, cfl, domlen, dt, time, timeout, timeto;
  int cells, N_max;
  vector<double> d, u, p, du, E, flux_d, flux_du, flux_E;

  // Initial conditions
  void computePrimitiveArraysToConservative();
  bool computeDirichletBoundaryConditions(const double& dl, 
					  const double& dul, 
					  const double& El, 
					  const double& dr, 
					  const double& dur, 
					  const double& Er);
  void computePeriodicBoundaryConditions();

  // Algorithm update
  void computeCFLCondition();
  void updateTimeStep();
  void computeIntercellFluxes();
  void updateCells();
  void outputCellData();

  // Output
  void computeConservativeArraysToPrimitive();
  void outputCellData();  

 public:
  CGodunov1DSolverEulerBase(const int& _cells, 
			    const int& _N_max, 
			    const double& _cfl, 
			    const double& _domlen, 
			    const double& _timeout,
			    const double& _gamma,
			    vector<double>& d, 
			    vector<double>& u, 
			    vector<double>& p);
  ~CGodunov1DSolverEulerBase();

  void computeSolution();
};

#endif
