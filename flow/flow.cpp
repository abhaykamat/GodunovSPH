#ifndef _FLOW_CPP
#define _FLOW_CPP


namespace flow {

  // StateVector
  // ===========

  StateVector::StateVector(const double& density,
			   const double& u_velocity,
			   const double& pressure) {
    state[0] = density;
    state[1] = density * u_velocity;  // U-Momentum
    state[2] = calcTotalPressureFromPrimitives(density, u_velocity, pressure); // Total Energy
  } 
  
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

  void StateVector::setPrimitives(const double& density, const double& u_velocity, const double& pressure) {
    state[0] = density;
    state[1] = density * u_velocity;
    state[2] = calcTotalPressureFromPrimitives(density, u_velocity, pressure);
  }

  void StateVector::setConservatives(const double& density, const double& u_momentum, const double& total_energy) {
    state[0] = density;
    state[1] = u_momentum;
    state[2] = total_energy;
  }

  // GammaConstants
  // ==============
  
  GammaConstants::GammaConstants(const double& _gamma) : gamma(_gamma) {
    g1 = (gamma - 1.0)/(2.0*gamma);
    g2 = (gamma + 1.0)/(2.0*gamma);
    g3 = 2.0*gamma/(gamma - 1.0);
    g4 = 2.0/(gamma - 1.0);
    g5 = 2.0/(gamma + 1.0);
    g6 = (gamma - 1.0)/(gamma + 1.0);
    g7 = (gamma - 1.0)/2.0;
    g8 = gamma - 1.0;
  }

  const double GammaConstants::getGamma() const { return gamma; }
  const double GammaConstants::getG1() const { return g1; }
  const double GammaConstants::getG2() const { return g2; }
  const double GammaConstants::getG3() const { return g3; }
  const double GammaConstants::getG4() const { return g4; }
  const double GammaConstants::getG5() const { return g5; }
  const double GammaConstants::getG6() const { return g6; }
  const double GammaConstants::getG7() const { return g7; }
  const double GammaConstants::getG8() const { return g8; }
  
}

// Compute the total energy, E, from primitive variables and ideal gas equation of state
double calcTotalEnergyFromPrimitives(const double& density, const double& u_velocity, const double& pressure) {
  // See p88 of [1] (Toro, 2nd Edition).
  return density*(0.5*u_velocity*u_velocity + pressure/((gamma - 1)*density));
}

// Calculate the primitive pressure variable from the conservative variables and ideal gas equation of state
void calcPressureFromConservatives(const double& density, const double& u_momentum, const double& total_energy) {
  return (gamma - 1)*(total_energy - 0.5*u_momentum*u_momentum/density);
}

#endif
