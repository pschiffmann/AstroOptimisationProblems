#include "trajobjfuns.h"
#include <iostream>
#include <vector>

int main() {
  // Values taken from the _solutions_ section of this website:
  // http://www.esa.int/gsp/ACT/inf/projects/gtop/gtoc1.html
  std::vector<double> X{6810.40521106, 168.37469758,  1079.47409963,
                        56.38731208,   1044.09288643, 3820.84181773,
                        1044.32726019, 3397.21349495};
  std::vector<double> rp;
  double objective_value = gtoc1(X, rp);

  for (double v : rp) {
    std::cout << v << ", " << std::endl;
  }
  std::cout << "Function Value: " << objective_value << std::endl;
}
