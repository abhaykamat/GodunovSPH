#ifndef _RIEMANN_CPP
#define _RIEMANN_CPP

#include "riemann.h"

namespace riemann {

  // RiemannSolverEulerBase
  // ======================

  RiemannSolverEulerBase::RiemannSolverEulerBase(const flow::StateVector& _left_state,
						 const flow::StateVector& _right_state) :
    left_state(_left_state), right_state(_right_state) {}

  RiemannSolverEulerBase::~RiemannSolverEulerBase() {}

  const flow::StateVector RiemannSolverEulerBase::getLeftState() const {
    return left_state;
  }

  const flow::StateVector RiemannSolverEulerBase::getRightState() const {
    return right_state;
  }
  
  // RiemannSolverEulerExact
  // =======================

  RiemannSolverEulerExact::RiemannSolverEulerExact(const flow::StateVector& _left_state, 
						   const flow::StateVector& _right_state) :
    RiemannSolverEulerBase(_left_state, _right_state) {
    init(_left_state, _right_state);
  }

  void RiemannSolverEulerExact::init(const flow::StateVector& _left_state,
				     const flow::StateVector& _right_state) {
    // Set convenience members for 
    // flow states on left and right
    dl = _left_state.getDensity();
    ul = _left_state.getUVelocity();
    pl = _left_state.getPressure();
    dr = _right_state.getDensity();
    ur = _right_state.getUVelocity();
    pr = _right_state.getPressure();

    // Initialise the Riemann Solver flow properties
    computeSoundSpeeds();
    bool test = testForVacuum();
    computePressureVelocityStar();
  }
  
  RiemannSolverEulerExact::~RiemannSolverEulerExact() {}

  void RiemannSolverEulerExact::computeSoundSpeeds() {
    cl = sqrt(flow::gammaconst::gamma * pl / dl);
    cr = sqrt(flow::gammaconst::gamma * pr / dr);
  }
  
  bool RiemannSolverEulerExact::testForVacuum() {
    if (flow::gammaconst::g4 * (cl + cr) <= (ur - ul)) {
      return true;
    }
    return false;
  }
  
  void RiemannSolverEulerExact::computePressureFunction(double& f, 
							double& fd, 
							const double& p, 
							const double& dk, 
							const double& pk, 
							const double& ck) {
    if (p <= pk) {   // Rarefaction wave
      double prat = p/pk;
      f = flow::gammaconst::g4 * ck * (pow(prat, flow::gammaconst::g1) - 1.0);
      fd = (1.0 / (dk * ck)) * pow(prat, -1.0*flow::gammaconst::g2);
    } else {        // Shock wave
      double ak = flow::gammaconst::g5 / dk;
      double bk = flow::gammaconst::g6 * pk;
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
    ppv = std::max(0.0, ppv);
    double p_min = std::min(pl, pr);
    double p_max = std::max(pl, pr);
    double q_max = p_max/p_min;
    
    if ((q_max <= q_user) && (p_min <= ppv && ppv <= p_max)) {   
      // Select PVRS Riemann solver
      p_iter = ppv;
    } else {
      if (ppv <= p_min) {   
	// Select Two-Rarefaction Riemann solver
	double pq = pow((pl / pr), flow::gammaconst::g1);
	double um = (pq*ul/cl + ur/cr + flow::gammaconst::g4*(pq - 1.0))/(pq/cl + 1.0/cr);
	double ptl = 1.0 + flow::gammaconst::g7*(ul - um)/cl;
	double ptr = 1.0 + flow::gammaconst::g7*(um - ur)/cr;
	p_iter = 0.5*(pl*pow(ptl, flow::gammaconst::g3) + 
		      pr*pow(ptr, flow::gammaconst::g3));
      } else {              
	// Select Two-Shock Riemann solver with PVRS as estimate
	double gel = sqrt((flow::gammaconst::g5/dl)/(flow::gammaconst::g6*pl + ppv));
	double ger = sqrt((flow::gammaconst::g5/dr)/(flow::gammaconst::g6*pr + ppv));
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
    double d, u, p;

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
	  double cml = (cl*pow((pm/pl),flow::gammaconst::g1));
	  double stl = um - cml;
	  if (wave_speed > stl) {
	    // Sampled point is Star Left state
	    d = dl*pow((pm/pl),(1.0/flow::gammaconst::gamma));
	    u = um;
	    p = pm;
	  } else {
	    // Sampled point is inside left fan
	    u = flow::gammaconst::g5*(cl + flow::gammaconst::g7*ul + wave_speed);
	    double c = flow::gammaconst::g5*(cl + flow::gammaconst::g7*(ul - wave_speed));
	    d = dl*pow((c/cl),flow::gammaconst::g4);
	    p = pl*pow((c/cl),flow::gammaconst::g3);
	  }
	}
      } else {
	// Left shock
	double pml = pm/pl;
	double sl = ul - cl*sqrt(flow::gammaconst::g2*pml + flow::gammaconst::g1);
	if (wave_speed <= sl) {
	  // Sampled point is left data state
	  d = dl;
	  u = ul;
	  p = pl;
	} else {
	  // Sampled point is Star Left state
	  d = dl * (pml + flow::gammaconst::g6)/(pml*flow::gammaconst::g6 + 1.0);
	  u = um;
	  p = pm;
	}
      }
    } else {
      // Sampling point lies to right of contact discontinuity
      if (pm > pr) {
	// Right shock
	double pmr = pm/pr;
	double sr = ur + cr*sqrt(flow::gammaconst::g2*pmr + flow::gammaconst::g1);
	if (wave_speed >= sr) {
	  // Sampled point is right data state
	  d = dr;
	  u = ur;
	  p = pr;
	} else {
	  // Sampled point is Star Right state
	  d = dr*(pmr + flow::gammaconst::g6)/(pmr*flow::gammaconst::g6 + 1.0);
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
	  double cmr = cr*pow((pm/pr),flow::gammaconst::g1);
	  double str = um + cmr;
	  if (wave_speed <= str) {
	    // Sampled point is Star Right state
	    d = dr*pow((pm/pr),(1.0/flow::gammaconst::gamma));
	    u = um;
	    p = pm;
	  } else {
	    u = flow::gammaconst::g5*(-1.0*cr + flow::gammaconst::g7*ur + wave_speed);
	    double c = flow::gammaconst::g5*(cr - flow::gammaconst::g7*(ur - wave_speed));
	    d = dr*pow((c/cr),flow::gammaconst::g4);
	    p = pr*pow((c/cr),flow::gammaconst::g3);
	  }
	}
      }
    }
    
    flow::StateVector state(d, u, p);
    return state;
  }

  void RiemannSolverEulerExact::testSolver(const double& domlen,
					   const double& diaph,
					   const int& cells,
					   const double& timeout) {
    std::ofstream outfile;
    outfile.open("/vol/godunovsph/test/riemann_test.out");

    double dx = domlen/static_cast<double>(cells);
    double x_pos, wave_speed;

    for (double i = 0; i < cells; i+=1.0) {
      x_pos = ((i+1.0) - 0.5)*dx;
      wave_speed = (x_pos - diaph)/timeout;
      
      flow::StateVector sol = sampleWaveSolution(wave_speed);
      outfile << x_pos << "\t" 
	      << sol.getDensity() << "\t" 
	      << sol.getUVelocity() << "\t" 
	      << sol.getPressure() << std::endl;
    }
    
    outfile.close();
  }
}

#endif
