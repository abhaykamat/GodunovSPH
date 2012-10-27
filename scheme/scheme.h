#ifndef _SCHEME_H
#define _SCHEME_H

#include "../mesh/mesh.h"
#include "../riemann/riemann.h"

namespace scheme {

  class ConservativeGRPScheme {
  private:
    double cfl, dt, time, time_out, time_to;
    mesh::Mesh1D flow_mesh;
    std::vector<flow::StateVector> fluxes;
    
    void setBoundaryConditions();  // Periodic for the time being
    void calcCFLCondition();
    void updateTimeStep();
    void calcIntercellFluxes();
    void updateCells();
    
  public:
    ConservativeGRPScheme(const double& _time_out, const int& _N, const double& _width);
    
    int setInitialConditions(const double& density_left, 
			     const double& momentum_left, 
			     const double& energy_left, 
			     const double& density_right, 
			     const double& momentum_right, 
			     const double& energy_right);
    void computeSolution();
    void outputSolution();
  };

}

#endif
