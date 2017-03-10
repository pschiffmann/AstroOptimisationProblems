#include "helpers.hpp"

std::array<double, 3> cross_product(std::array<double, 3> v1,
                                    std::array<double, 3> v2) {
  return {v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2],
          v1[0] * v2[1] - v1[1] * v2[0]};
}

double dot_product(std::array<double, 3> v1, std::array<double, 3> v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
