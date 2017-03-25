#include "vector3d_helpers.hpp"
#include <cmath>

std::array<double, 3> cross_product(std::array<double, 3> v1,
                                    std::array<double, 3> v2) {
  return {v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2],
          v1[0] * v2[1] - v1[1] * v2[0]};
}

double dot_product(std::array<double, 3> v1, std::array<double, 3> v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

double norm(std::array<double, 3> v) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

std::array<double, 3> unit_vector(const std::array<double, 3>& v) {
  return div(v, norm(v));
}

std::array<double, 3> add(std::array<double, 3> v1, std::array<double, 3> v2) {
  return {v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]};
}

std::array<double, 3> sub(std::array<double, 3> v1, std::array<double, 3> v2) {
  return {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
}

std::array<double, 3> mul(std::array<double, 3> v1, const double r) {
  return {v1[0] * r, v1[1] * r, v1[2] * r};
}

std::array<double, 3> div(std::array<double, 3> v1, const double r) {
  return {v1[0] / r, v1[1] / r, v1[2] / r};
}
