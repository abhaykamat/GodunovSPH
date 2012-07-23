#ifndef _RIEMANN_CPP
#define _RIEMANN_CPP

#include "riemann.h"

namespace riemann {

  // RiemannSolverEulerBase
  // ======================

  RiemannSolverEulerBase::RiemannSolverEulerBase(const flow::StateVector& _left_state,
						 const flow::StateVector& _right_state) :
    left_state(_left_state), right_state(_right_state) {}

  const flow::StateVector RiemannSolverEulerBase::getLeftSate() {
    return left_state;
  }

  const flow::StateVector RiemannSolverEulerBase::getRightState() {
    return right_state;
  }
  
  // RiemannSolverEulerExact
  // =======================
  
  RiemannSolverEulerExact::RiemannSolverEulerExact(const flow::StateVector& _left_state, const flow::StateVector& _right_state) : RiemannSolverEulerBase(_left_state, _right_state) {}
  
  void RiemannSolverEulerExact::computeSoundSpeeds() {
    left_speed_of_sound = sqrt(gamma * pl / dl);
    right_speed_of_sound = sqrt(gamma * pr / dr);
  }
  
  bool RiemannSolverEulerExact::testForVacuum() {
    if (g4 * (cl + cr) <= (ur - ul)) {
      return true;
    }
    return false;
  }
  
  void RiemannSolverEulerExact::computePressureFunction(double& f, double& fd, const double& p, const double& dk, const double& pk, const double& ck) {
    if (p <= pk) {   // Rarefaction wave
      double prat = p/pk;
      f = g4 * ck * (pow(prat, g1) - 1.0);
      fd = (1.0 / (dk * ck)) * pow(prat, -1.0*g2);
    } else {        // Shock wave
      double ak = g5 / dk;
      double bk = g6 * pk;
      double qrt = sqrt( ak / (bk + p));
      f = (p - pk) * qrt;
      fd = (1.0 - 0.5 * (p - pk) / (bk + p)) * qrt;
    }
  }
  
  double RiemannSolverEulerExact::computeGuessPressure(const double& p_start) {
    double p_iter = p_start;
    double q_user = 2.0;
    
    // Compute guess pressure from PVRS Riemann Solver
    double cup = 0.25 * (dl + dr) * (cl + cr);
    double ppv = 0.5 * (pl + pr) + 0.5 * (ul - ur) * cup;
    ppv = max(0.0, ppv);
    double p_min = min(pl, pr);
    double p_max = max(pl, pr);
    double q_max = p_max/p_min;
    
    if ((q_max <= q_user) && (p_min <= ppv && ppv <= p_max)) {   
      // Select PVRS Riemann solver
      p_iter = ppv;
    } else {
      if (ppv <= p_min) {   
	// Select Two-Rarefaction Riemann solver
	double pq = pow((pl / pr), g1);
	double um = (pq*ul/cl + ur/cr + g4*(pq - 1.0))/(pq/cl + 1.0/cr);
	double ptl = 1.0 + g7*(ul - um)/cl;
	double ptr = 1.0 + g7*(um - ur)/cr;
	p_iter = 0.5*(pl*pow(ptl, g3) + pr*pow(ptr, g3));
      } else {              
	// Select Two-Shock Riemann solver with PVRS as estimate
	double gel = sqrt((g5/dl)/(g6*pl + ppv));
	double ger = sqrt((g5/dr)/(g6*pr + ppv));
	p_iter = (gel*pl + ger*pr - (ur - ul))/(gel + ger); 
      }
    }
    
    return p_iter;
  }
  
  void RiemannSolverEulerExact::computePressureVelocityStar() {
    const double tol_pre = 1.0e-6;
    const int nr_iter = 20;
    
    const double p_start = 0.0;
    double p_old = computeGuessPressure(p_start);
    
    const double u_diff = ur - ul;
    double fl, fr, fld, frd = 0.0;
    double change = 0.0;
    
    // Compute pressure in star region
    for (int i = 0; i < nr_iter; i++) {
      computePressureFunction(fl, fld, p_old, dl, pl, cl);
      computePressureFunction(fr, frd, p_old, dr, pr, cr);
      pm = p_old - (fl + fr + u_diff)/(fld + frd);
      change = 2.0 * fabs((pm - p_old)/(pm + p_old));
      if (change <= tol_pre) {
	break;
      }
      if (pm < 0.0) {
	pm = tol_pre;
      }
      p_old = pm;
    }
    
    // Compute velocity in star region
    um = 0.5*(ul + ur + fr - fl);
  }
  
  flow::StateVector RiemannSolverEulerExact::sampleWaveSolution(const double& wave_speed) {
    if (wave_speed <= um) {
      // Sampling point lies to the left of the contact discontinuity
      if (pm <= pl) {
	// Left rarefaction
	double shl = ul - cl;
	if (wave_speed <= shl) {
	  // Sampled point is left data state
	  d = dl;
	  u = ul;
	  p = pl;
	} else {
	  double cml = (cl*pow((pm/pl),g1));
	  double stl = um - cml;
	  if (wave_speed > stl) {
	    // Sampled point is Star Left state
	    d = dl*pow((pm/pl),(1.0/gamma));
	    u = um;
	    p = pm;
	  } else {
	    // Sampled point is inside left fan
	    u = g5*(cl + g7*ul + s);
	    double c = g5*(cl + g7*(ul - s));
	    d = dl*pow((c/cl),g4);
	    p = pl*pow((c/cl),g3);
	  }
	}
      } else {
	// Left shock
	double pml = pm/pl;
	double sl = ul - cl*sqrt(g2*pml + g1);
	if (wave_speed <= sl) {
	  // Sampled point is left data state
	  d = dl;
	  u = ul;
	  p = pl;
	} else {
	  // Sampled point is Star Left state
	  d = dl * (pml + g6)/(pml*g6 + 1.0);
	  u = um;
	  p = pm;
	}
      }
    } else {
      // Sampling point lies to right of contact discontinuity
      if (pm > pr) {
	// Right shock
	double pmr = pm/pr;
	double sr = ur + cr*sqrt(g2*pmr + g1);
	if (wave_speed >= sr) {
	  // Sampled point is right data state
	  d = dr;
	  u = ur;
	  p = pr;
	} else {
	  // Sampled point is Star Right state
	  d = dr*(pmr + g6)/(pmr*g6 + 1.0);
	  u = um;
	  p = pm;
	}
      } else {
	// Right rarefaction
	double shr = ur + cr;
	if (wave_speed >= shr) {
	  // Sampled point is right data state
	  d = dr;
	  u = ur;
	  p = pr;
	} else {
	  double cmr = cr*pow((pm/pr),g1);
	  double str = um + cmr;
	  if (s <= str) {
	    // Sampled point is Star Right state
	    d = dr*pow((pm/pr),(1.0/gamma));
	    u = um;
	    p = pm;
	  } else {
	    u = g5*(-1.0*cr + g7*ur + s);
	    double c = g5*(cr - g7*(ur - s));
	    d = dr*pow((c/cr),g4);
	    p = pr*pow((c/cr),g3);
	  }
	}
      }
    }
    
    flow::StateVector state(d, u, p);
    return state;
  }
  
}

#endif
