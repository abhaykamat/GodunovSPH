#include "scheme.h"

int main() {
  ConservativeGRPScheme scheme = ConservativeGRPScheme(0.25, 100, 1.0);
  scheme.computeSolution();
  scheme.outputSolution();
  return 0;
}
