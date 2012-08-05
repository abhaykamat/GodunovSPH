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

  StateVector ConservativeGRPScheme::calcAcousticDerivative(const RiemannSolverEulerExact& rs, const StateVector& U_jphalf, const int& j) {
    // See pg 178-179 of [2] for full mathematical treatment

    // (A) If r = r_{j+1/2} is to the right of \Gamma_3 or to the left
    // of \Gamma_3, acoustic derivative is determined by initial data
    bool lor = rs.leftOrRightWaveTrue();
    if (lor == true) {
      StateVector acoustic = euler::calcEulerFluxFromConservatives(U_jphalf.getDensity(),
								   U_jphalf.getUMomentum(),
								   U_jphalf.getPressure());
      return StateVector acous;
    }

    // Suppose that r = r_{j+1/2} is between \Gamma_1 and \Gamma_3

    // Redefinition of left states for ease of following algorithm
    double rhoL = mesh.meshPoint(j).getDensity();
    double uL = mesh.meshPoint(j).getUVelocity();
    double pL = mesh.meshPoint(j).getPressure();
    double cL = sqrt(flow::gammaconst::gamma * pL / rhoL);
    double gL = rhoL*cL;

    double drhodxL = mesh.slopePoint(j).getDensity();
    double dudxL = mesh.slopePoint(j).getUVelocity();
    double dpdxL = mesh.slopePoint(j).getPressure();
    
    // Redefinition of right states for ease of following algorithm
    double rhoR = mesh.meshPoint(j+1).getDensity();
    double uR = mesh.meshPoint(j+1).getUVelocity();
    double pR = mesh.meshPoint(j+1).getPressure();
    double cR = sqrt(flow::gammaconst::gamma * pR / rhoR);
    double gR = rhoR*cR;

    double dudxR = mesh.slopePoint(j+1).getDensity();
    double dudxR = mesh.slopePoint(j+1).getUVelocity();
    double dpdxR = mesh.slopePoint(j+1).getPressure();

    // TODO: Need to obtain these values (rhoQStar, cQStar)!
    double lambda = 1.0;
    double gLStar = rhoLStar * cLStar; 
    double gRStar = rhoRStar * cRStar;

    // A, B and D coefficients to find star time derivatives
    aL = 1.0;
    aR = -1.0;

    bL = 1.0/gLStar;
    bR = 1.0/gRStar;

    dL = -1.0*bL*(cL*gL*dudxL + cL*dpdxL - lambda*uL*cL*gL);
    dR = -1.0*bR*(cR*gR*dudxR - cR*dpdxR + lambda*uR*cR*gR);

    // Solve linear simultaneous equations to find time derivates 
    // TODO: Fill in the values for dudtStar, dpdtStar
    double dudtStar = 0.0;
    double dpdtStar = 0.0;

    double drhodtLStar = (1.0/(cLStar*cLStar)) * dpdtStar;
    double drhodtRStar = (1.0/(cRStar*cRStar)) * dpdtStar;

    // Section (D): Evaluate the derivatives at the contact discontinuity
    double dpdxiStar = -1.0*dudtStar;
    double dudxiRStar = -1.0*(1.0/(gRStar*gRStar))*dpdtStar - lambda*(1.0/rhoRStar)*uStar;
    double dudxiLStar = -1.0*(1.0/(gLStar*gLStar))*dpdtStar - lambda*(1.0/rhoLStar)*uStar;
    double drhodxiLStar = (1.0/gLStar)*(rhoL*dudxL + cL*drhodxL + lambda*rhoL*uL + 
					(1.0/(cLStar*cLStar))*dpdtStar);
    double drhodxiRStar = (1.0/gRStar)*(-1.0*rhoR*dudxR + cR*drhodxR - lambda*rhoR*uR -
					(1.0/(cRStar*cRStar))*dpdtStar);

    // Section (E): Evaulate the final acoustic time derivative
    // TODO: Determine whether sample is between \Gamma1 & \Gamma2 or \Gamma2 and \Gamma3
    // if between gamma1 and gamma2
    StateVector dUdt0(drhodtLStar, dudtStar, dpdtStar);
    StateVector dUdxi0(drhodxiLStar, dudxiLStar, dpdxiStar);
    // if between gamma2 and gamma3
    StateVector dUdt0(drhodtRStar, dudtStar, dpdtStar);
    StateVector dUdxi0(drhodxiRStar, dudxiRStar, dpdxiStar);
    
    // TODO: Will require operator* for scalar values, operator+ for StateVector values
    // TODO: Need to obtain rho0 and u0 (these are from the U_jphalf)
    StateVector acoustic = dUdt0 - rho0*u0*dUdxi0;
    return acoustic;
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
      
      // Calculate the acoustic derivatives 
      // (\frac{\partial}{\partial t} U^{ac})^n_{j+1/2}), as for pg 178-179 of [2]
      StateVector ddt_U_ac_jphalf = calcAcousticDerivative(rs, U_jphalf, i);

      // Update the fluxes. This is the E_1 scheme.
      // TODO: Determine the value of k.
      double k_half = k/2.0;
      double flux_d = U_jphalf.getDensity() + k_half * ddt_U_ac_jphalf.getDensity();
      double flux_um = U_jphalf.getUMomentum() + k_half * ddt_U_ac_jphalf.getUMomentum();
      double flux_p = U_jphalf.getTotalEnergy() + k_half *ddt_U_ac_jphalf.getTotalEnergy();

      fluxes.meshPoint(i) = euler::calcEulerFluxFromConservatives(flux_d, flux_um, flux_p);
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

  ConservativeGRPScheme::ConservativeGRPScheme(const double& _time_out, 
					       const int& _N, 
					       const double& _width) : 
    time_out(_time_out), mesh(_N, _width) {
    
  }

  void ConservativeGRPScheme::computeSolution() {

  }

  void ConservativeGRPScheme::outputSolution() {

  }

}

#endif
