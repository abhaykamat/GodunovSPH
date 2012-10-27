#include "mesh/mesh.h"

int main() {
  //scheme::ConservativeGRPScheme grpscheme = scheme::ConservativeGRPScheme(0.25, 100, 1.0);
  //scheme.computeSolution();
  //scheme.outputSolution();
  return 0;
}

/*
#include "./riemann/riemann.h"

int main() {
  flow::StateVector left = flow::StateVector(1.0, 0.0, 1.0);
  flow::StateVector right = flow::StateVector(0.125, 0.0, 0.1);
  riemann::RiemannSolverEulerExact rs = riemann::RiemannSolverEulerExact(left, right);
  rs.testSolver(1.0, 0.5, 100, 0.25);
  return 0;
  }*/
