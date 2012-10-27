#ifndef _SCHEME_CPP
#define _SCHEME_CPP

#include "scheme.h"

namespace scheme {

  void ConservativeGRPScheme::setBoundaryConditions() {
    int N = flow_mesh.size();
    flow_mesh.meshPoint(N) = flow_mesh.meshPoint(1);
    flow_mesh.meshPoint(1) = flow_mesh.meshPoint(N-1);
  }

  int ConservativeGRPScheme::setInitialConditions(const double& density_left, 
						  const double& momentum_left, 
						  const double& energy_left, 
						  const double& density_right, 
						  const double& momentum_right, 
						  const double& energy_right) {
    if (density_left < 0.0 || energy_left < 0.0 || 
	density_right < 0 || energy_right < 0) {   // Unphysical conditions
      return 0;
    }
    
    int N_half = fabs(static_cast<double>(flow_mesh.size())/2.0);
    for (int i=1; i<N_half; i++) {
      flow_mesh.meshPoint(i).setPrimitives(density_left, momentum_left, energy_left);
    }
    for (int i=N_half; i<=flow_mesh.size(); i++) {
      flow_mesh.meshPoint(i).setPrimitives(density_right, momentum_right, energy_right);
    }

    return 1;
  }

  void ConservativeGRPScheme::calcCFLCondition() {
      double smax = -1.0e6;
      int N = flow_mesh.size();

      // Find the maximum characteristic speed
      for (int i = 1; i <= N; i++) {
	double du_d = fabs(flow_mesh.meshPoint(i).getUMomentum() / flow_mesh.meshPoint(i).getDensity());

	if (du_d > smax) {
	  smax = du_d;
	}
      }
      
      // Restrict the time step via CFL
      dt = cfl * flow_mesh.cell_size(1) / smax;
      
      // Check size of DT to avoid exceeding output time
      if ((time + dt) > time_out) {
	// Recompute DT
	dt = time_out - time;
      }
  }

  void ConservativeGRPScheme::updateTimeStep() {
    time += dt;
  }

  /*
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
    }*/

  void ConservativeGRPScheme::calcIntercellFluxes() {
    for (int i = 1; i <= flow_mesh.size(); i++) {
      // Calculate the updated boundary values as for pg 174 of [2]
      //StateVector dx_half = mesh.dx(i)/2.0 * mesh.slopePoint(i);
      flow::StateVector U_jphalf_minus = flow_mesh.meshPoint(i);// + dx_half;
      flow::StateVector U_jphalf_plus = flow_mesh.meshPoint(i);// - dx_half;

      // Use the Riemann Solver with the updated boundary values
      riemann::RiemannSolverEulerExact rs = riemann::RiemannSolverEulerExact(U_jphalf_minus, U_jphalf_plus);
      flow::StateVector U_jphalf = rs.sampleWaveSolution(0.0);
      
      // Calculate the acoustic derivatives 
      // (\frac{\partial}{\partial t} U^{ac})^n_{j+1/2}), as for pg 178-179 of [2]
      //StateVector ddt_U_ac_jphalf = calcAcousticDerivative(rs, U_jphalf, i);

      // Update the fluxes. This is the E_1 scheme.
      // TODO: Determine the value of k.
      //double k_half = k/2.0;
      double flux_d = U_jphalf.getDensity();// + k_half * ddt_U_ac_jphalf.getDensity();
      double flux_um = U_jphalf.getUMomentum();// + k_half * ddt_U_ac_jphalf.getUMomentum();
      double flux_E = U_jphalf.getTotalEnergy();// + k_half *ddt_U_ac_jphalf.getTotalEnergy();

      fluxes[i] = flow::calcEulerFluxFromConservatives(flux_d, flux_um, flux_E);
    }
  }

  void ConservativeGRPScheme::updateCells() {
    for (int i = 1; i <= flow_mesh.size(); i++) {
      double dtodx = dt/flow_mesh.cell_size(i);

      double density = flow_mesh.meshPoint(i).getDensity() + 
	dtodx * (fluxes[i-1].getDensity() - fluxes[i].getDensity());
      double u_momentum = flow_mesh.meshPoint(i).getUMomentum() + 
	dtodx * (fluxes[i-1].getUMomentum() - fluxes[i].getUMomentum());
      double total_energy = flow_mesh.meshPoint(i).getTotalEnergy() + 
	dtodx * (fluxes[i-1].getTotalEnergy() - fluxes[i].getTotalEnergy());

      flow_mesh.meshPoint(i).setConservatives(density, u_momentum, total_energy);
    }
  }

  ConservativeGRPScheme::ConservativeGRPScheme(const double& _time_out, 
					       const int& _N, 
					       const double& _width) : 
    time_out(_time_out), flow_mesh(_N, _width) {
  }

  void ConservativeGRPScheme::computeSolution() {
    setBoundaryConditions();     // Use periodic boundary conditions (for now)
    calcCFLCondition();       // Impose the Courant-Friedrichs-Lewy (CFL) condition
    updateTimeStep();            // Update the global time step
    calcIntercellFluxes();    // Compute intercell numerical fluxes
    updateCells();               // Update the solution according to conservative formula
  }

  void ConservativeGRPScheme::outputSolution() {
    std::ofstream outfile;
    outfile.open("/vol/godunovsph/test/scheme_test.out");
    double x_pos = 0.0;
    for (double i = 0; i < flow_mesh.size(); i+=1.0) {
      x_pos = ((i+1.0) - 0.5)*flow_mesh.cell_size(i);
      outfile << x_pos << "\t" 
	      << flow_mesh.meshPoint(i).getDensity() << "\t" 
	      << flow_mesh.meshPoint(i).getUVelocity() << "\t" 
	      << flow_mesh.meshPoint(i).getPressure() << std::endl;
    }
    
    outfile.close();
  }

}

#endif
