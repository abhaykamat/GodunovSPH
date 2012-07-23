#ifndef _SCHEME_CPP
#define _SCHEME_CPP

#include "scheme.h"

namespace scheme {

  void ConservativeGRPScheme::setInitialConditions() {

  }

  void ConservativeGRPScheme::calcCFLCondition() {

  }

  void ConservativeGRPScheme::updateTimeStep() {
    time += dt;
  }

  StateVector ConservativeGRPScheme::calcAcousticDerivative(const int& j) {
    // TODO: Lots of work in here!
  }

  void ConservativeGRPScheme::calcIntercellFluxes() {
    for (int i = 1; i <= N; i++) {

      // Calculate the updated boundary values as for pg 174 of [2]
      StateVector dx_half = mesh.dx(i)/2.0 * mesh.slopePoint(i);
      StateVector U_jphalf_minus = mesh.meshPoint(i) + dx_half;
      StateVector U_jphalf_plus = mesh.meshPoint(i) - dx_half;

      // Use the Riemann Solver with the updated boundary values
      RiemannSolverEulerExact rs = RiemannSolverEulerExact(U_jphalf_minus, U_jphalf_plus);
      StateVector U_jphalf = rs.sampleWaveSolution(0.0);
      
      // Calculate the acoustic derivatives (\frac{\partial}{\partial t} U^{ac})^n_{j+1/2}), as for pg 178-179 of [2]
      StateVector ddt_U_ac_jphalf = calcAcousticDerivative(i);

      // Update the fluxes (setting k = ...?). This is the E_1 scheme.
      double k_half = k/2.0;
      fluxes.meshPoint(i).setConservatives(
        U_jphalf.getDensity() + k_half * ddt_U_ac_jphalf.getDensity(),
	U_jphalf.getUMomentum() + k_half * ddt_U_ac_jphalf.getUMomentum(),
	U_jphalf.getTotalEnergy() + k_half *ddt_U_ac_jphalf.getTotalEnergy()
      );
    }
  }

  void ConservativeGRPScheme::updateCells() {
    for (int i = 1; i <= N; i++) {
      double density = mesh.meshPoint(i).getDensity() + 
	dtodx * (fluxes.meshPoint(i-1).getDensity() - 
		 fluxes.meshPoint(i).getDensity());
      double u_momentum = mesh.meshPoint(i).getUMomentum() + 
	dtodx * (fluxes.meshPoint(i-1).getUMomentum() - 
		 fluxes.meshPoint(i).getUMomentum());
      double total_energy = mesh.meshPoint(i).getTotalEnergy() + 
	dtodx * (fluxes.meshPoint(i-1).getTotalEnergy() - 
		 fluxes.meshPoint(i).getTotalEnergy());

      mesh.meshPoint(i).setConservatives(density, u_momentum, total_energy);
    }
  }

  ConservativeGRPScheme::ConservativeGRPScheme(const double& _time_out, const int& _N, const double& _width) : time_out(_time_out), mesh(_N, _width) {
    
  }

  void ConservativeGRPScheme::computeSolution() {

  }

  void ConservativeGRPScheme::outputSolution() {

  }

}

#endif
