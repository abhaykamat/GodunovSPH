#ifndef _SCHEME_H
#define _SCHEME_H

#include "riemann.h"

namespace scheme {

  class ConservativeGRPScheme {
    private:
      double cfl, dt, time, time_out, time_to;
      mesh::Mesh1D mesh;
      vector<flow::StateVector> fluxes;

      void setInitialConditions();
      void calcCFLCondition();
      void updateTimeStep();
      void calcIntercellFluxes();
      void updateCells();
      
  public:
      ConservativeGRPScheme(const double& _time_out, Mesh1D& _mesh);

      void computeSolution();
      void outputSolution();
  }

}

#endif
