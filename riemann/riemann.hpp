#ifndef _RIEMANN_HPP
#define _RIEMANN_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

class CRiemannSolverEulerBase {
protected:
  double dl, ul, pl, dr, ur, pr;
  double gamma, g1, g2, g3, g4, g5, g6, g7, g8;
  
public:
  CRiemannSolverEulerBase(const double& _dl,
			  const double& _ul,
			  const double& _pl,
			  const double& _dr,
			  const double& _ur,
			  const double& _pr,
			  const double& _gamma);
  ~CRiemannSolverEulerBase();
  
  void computeGammaConstants();
  
  virtual void sampleWaveSolution(const double& s,    // Wave speed, s = x/t 
				  double& d,          // Sampled density
				  double& u,          // Sampled velocity
				  double& p) = 0;     // Sampled pressure

  virtual void testSolver(const int& domlen,            // Domain length
			  const double& diaph,          // Initial discontinuity position
			  const int& cells,             // Number of computing cells
			  const double& timeout) = 0;   // Output time
			  };

class CRiemannSolverEulerExact: public CRiemannSolverEulerBase {
private:

  // Sound speeds
  double cl;
  double cr;

  // Pressure/Velocity in Star region
  double pm;
  double um;

  void computeSoundSpeeds();

  bool testForVacuum();

  void computePressureFunction(double& f,
			       double& fd,
			       const double& p,
			       const double& dk,
			       const double& pk,
			       const double& ck);

  void computeGuessPressure(double& p_start);

  void computePressureVelocityStar();

public:
  CRiemannSolverEulerExact(const double& _dl,    // Left density state
			   const double& _ul,    // Left velocity state
			   const double& _pl,    // Left pressure state
			   const double& _dr,    // Right density state
			   const double& _ur,    // Right velocity state
			   const double& _pr,    // Right pressure stat
			   const double& _gamma);   // Ratio of specific heats
  ~CRiemannSolverEulerExact();

  virtual void sampleWaveSolution(const double& s,
                                  double& d,
                                  double& u,
                                  double& p);

  virtual void testSolver(const int& domlen,     // Domain length
                          const double&diaph,    // Initial discontinuity position
                          const int& cells,      // Number of computing cells
                          const double&timeout); // Output time    
};

#endif
