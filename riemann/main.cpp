#include "riemann.hpp"

int main() {
  CRiemannSolverEulerExact rs = CRiemannSolverEulerExact(1.0, 0.0, 1.0, 0.125, 0.0, 0.1, 1.4);
  rs.testSolver(1.0, 0.5, 100, 0.25);
  return 0;
}
