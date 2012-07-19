#ifndef _RIEMANN_H
#define _RIEMANN_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>


class StateVector {
 private:
  std::vector<double> state(3);

 public:
  explicit StateVector(const double& density,
		       const double& velocity,
		       const double& pressure);
  virtual ~StateVector();

  const double getDensity() const;
  const double getVelocity() const;
  const double getPressure() const; 
};


class RiemannSolverEulerBase {
 private:
  StateVector left_state;
  StateVector right_state;

 public:
  explicit RiemannSolverEulerBase(const StateVector& _left_state,
				  const StateVector& _right_state);
  virtual ~RiemannSolverEulerBase();

  const StateVector getLeftState() const;
  const StateVector getRightState() const;
  
  virtual StateVector sampleWaveSolution(const double& wave_speed) = 0;
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

  virtual StateVector sampleWaveSolution(const double& wave_speed);
};

#endif
