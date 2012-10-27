#ifndef _MESH_H
#define _MESH_H

#include "../flow/flow.h"

namespace mesh {
  
  class Mesh1D {
    private:
      int N;                             // Number of cells
      double width;                      // Width of domain
      std::vector<double> dx;            // Width of cells
      std::vector<flow::StateVector> U;  // Flow value vector (conservative formualation)
      std::vector<flow::StateVector> L;  // Slope value vector (conservative formulation)
      
      void init();                       // Initialise the cell size vector

    public:
      Mesh1D(int _N, const double& _width);

      int size() const;
      double cell_size(const int& i) const;
      flow::StateVector& meshPoint(const int& j);
      flow::StateVector& slopePoint(const int& j);
  };

}

#endif
