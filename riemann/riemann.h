#ifndef _RIEMANN_H
#define _RIEMANN_H

#include <cmath>
#include <algorithm>
#include <fstream>
#include "../flow/flow.h"

namespace riemann {

  class RiemannSolverEulerBase {
    private:
      flow::StateVector left_state;
      flow::StateVector right_state;

    public:
      RiemannSolverEulerBase(const flow::StateVector& _left_state,
			     const flow::StateVector& _right_state);
      virtual ~RiemannSolverEulerBase();

      const flow::StateVector getLeftState() const;
      const flow::StateVector getRightState() const;
  
      virtual flow::StateVector sampleWaveSolution(const double& wave_speed) = 0;
  };

  class RiemannSolverEulerExact : public RiemannSolverEulerBase {
    private:
      // Convenience members for left and right states 
      double dl, ul, pl, dr, ur, pr;

      // Sounds speeds
      double cl, cr;

      // Pressure and Velocity in the Star Region
      double pm, um;

      void init(const flow::StateVector& _left_state, 
		const flow::StateVector& _right_state);
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
      RiemannSolverEulerExact(const flow::StateVector& _left_state,
			      const flow::StateVector& _right_state);
      virtual ~RiemannSolverEulerExact();

      virtual flow::StateVector sampleWaveSolution(const double& wave_speed);
      void testSolver(const double& domlen,
		      const double& diaph,
		      const int& cells,
		      const double& timeout);
  };

}
#endif
