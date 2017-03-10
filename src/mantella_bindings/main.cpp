#include <array>
#include <iostream>
#include "gtoc1.hpp"

int main() {
  multiple_gravity_assist::gtoc1 problem;

  // Values taken from the _solutions_ section of this website:
  // http://www.esa.int/gsp/ACT/inf/projects/gtop/gtoc1.html
  std::array<double, 8> flyby_sequence{
      6810.40521106, 168.37469758,  1079.47409963, 56.38731208,
      1044.09288643, 3820.84181773, 1044.32726019, 3397.21349495};
  double objective_value = problem.objective_function(flyby_sequence);

  std::cout << "Function Value: " << objective_value << std::endl;
}
