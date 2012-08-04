#ifndef _MESH_CPP
#define _MESH_CPP

#include "mesh.h"

namespace mesh {

  Mesh1D::init() {
    double cell_size = width/static_cast<double>(N);

    // Initialise the cell sizes with a constant value, cell_size
    for(std::vector<StateVector>::iterator it = dx.begin(); it != dx.end(); ++it) {
      *it = cell_size;
    }
  }

  Mesh1D::Mesh1D(size_t _N, const double& _width) : N(_N), width(_width), dx(_N, 0.0), U(_N, 0.0), L(_N, 0.0) {
    init();
  }

  StateVector& Mesh1D::meshPoint(const int& j) {
    return U[j-1];
  }

  StateVector& Mesh1D::slopePoint(const int& j) {
    return L[j-1];
  }

}

#endif


