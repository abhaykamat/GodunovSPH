#ifndef _FLOW_CPP
#define _FLOW_CPP

#include "riemann.h"

namespace flow {

  // StateVector
  // ===========

  StateVector::StateVector(const double& density,
			   const double& velocity,
			   const double& pressure) {
    state[0] = density;
    state[1] = velocity;
    state[2] = pressure;
  } 
  
  const double StateVector::getDensity() {
    return state[0];
  }
  
  const double StateVector::getVelocity() {
    return state[1];
  }
  
  const double StateVector::getPressure() {
    return state[2];
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

#endif
