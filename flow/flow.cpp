#ifndef _FLOW_CPP
#define _FLOW_CPP

#include "flow.h"

namespace flow {

  // Pressure/Energy Conversions
  // ===========================

  // Compute the total energy, E, from primitive variables 
  // and ideal gas equation of state
  double calcTotalEnergyFromPrimitives(const double& density, 
				       const double& u_velocity, 
				       const double& pressure) {
    // See p88 of [1] (Toro, 2nd Edition).
    return density*(0.5*u_velocity*u_velocity + 
		    pressure/((flow::gammaconst::gamma - 1)*density));
  }

  // Calculate the primitive pressure variable from 
  // the conservative variables and ideal gas equation of state
  double calcPressureFromConservatives(const double& density, 
				       const double& u_momentum, 
				       const double& total_energy) {
    return (flow::gammaconst::gamma - 1) * 
      (total_energy - 0.5*u_momentum*u_momentum/density);
  }

  // StateVector
  // ===========

  StateVector::StateVector(const double& density,
			   const double& u_velocity,
			   const double& pressure) {
    state.resize(3);
    state[0] = density;
    state[1] = density * u_velocity;  // U-Momentum
    state[2] = calcTotalEnergyFromPrimitives(density, u_velocity, pressure); // Total Energy
  } 

  StateVector::~StateVector() {}
  
  const double StateVector::getDensity() const {
    return state[0];
  }
  
  const double StateVector::getUVelocity() const {
    return state[1]/state[0];
  }
  
  const double StateVector::getUMomentum() const {
    return state[1];
  }

  const double StateVector::getPressure() const {
    return calcPressureFromConservatives(state[0], state[1], state[2]);
  }

  const double StateVector::getTotalEnergy() const {
    return state[2];
  }

  void StateVector::setPrimitives(const double& density, 
				  const double& u_velocity, 
				  const double& pressure) {
    state[0] = density;
    state[1] = density * u_velocity;
    state[2] = calcTotalEnergyFromPrimitives(density, u_velocity, pressure);
  }

  void StateVector::setConservatives(const double& density, 
				     const double& u_momentum, 
				     const double& total_energy) {
    state[0] = density;
    state[1] = u_momentum;
    state[2] = total_energy;
  }  

}

#endif
