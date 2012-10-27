#include "scheme/scheme.h"
#include <iostream>

int main() {
  scheme::ConservativeGRPScheme grpscheme = scheme::ConservativeGRPScheme(0.25, 100, 1.0);
  scheme.setBoundaryConditions();
  scheme.computeSolution();
  scheme.outputSolution();
  //mesh::Mesh1D flow_mesh(100, 1.0);
  //std::cout << flow_mesh.cell_size(1) << std::endl;
  return 0;
}
