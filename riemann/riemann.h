#ifndef _RIEMANN_H
#define _RIEMANN_H

#include <cmath>
#include <algorithm>
#include "flow/flow.h"


class RiemannSolverEulerBase {
 private:
  flow::StateVector left_state;
  flow::StateVector right_state;

 public:
  explicit RiemannSolverEulerBase(const StateVector& _left_state,
				  const StateVector& _right_state);
  virtual ~RiemannSolverEulerBase();

  const flow::StateVector getLeftState() const;
  const flow::StateVector getRightState() const;
  
  virtual flow::StateVector sampleWaveSolution(const double& wave_speed) = 0;
};


class RiemannSolverEulerExact : public RiemannSolverEulerBase {
 private:

  // Sounds speeds
  double left_speed_of_sound;
  double right_speed_of_sound;

  // Pressure and Velocity in the Star Region
  double pressure_star;
  double velocity_star;

  void computeSoundSpeeds();
  bool testForVacuum();
  void computePressureFunction(double& f,
			       double& fd,
			       const double& p,
			       const double& dk,
			       const double& pk,
			       const double& ck);
  double computeGuessPressure(const double& p_start);
  void computePressureVelocityStar();
  
 public:
  explicit RiemannSolverEulerExact();
  virtual ~RiemannSolverEulerExact();

  virtual flow::StateVector sampleWaveSolution(const double& wave_speed);
};

#endif
