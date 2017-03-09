#pragma once

#include <array>

#include <mantella0>

namespace gtoc1 {
struct celestial_body {
 public:
  static const celestial_body SUN;
  static const celestial_body MERCURY;
  static const celestial_body VENUS;
  static const celestial_body EARTH;
  static const celestial_body MARS;
  static const celestial_body JUPITER;
  static const celestial_body SATURN;
  static const celestial_body URANUS;
  static const celestial_body NEPTUNE;

 private:
  celestial_body(double mu, double penalty,
                 double penalty_coefficient) noexcept;
  double mu;
  double penalty;
  double penalty_coefficient;
};

struct asteroid {
  asteroid(std::array<double, 6> keplerian, double epoch, double mu) noexcept;

  std::array<double, 6> keplerian;
  double epoch;
  double mu;
};

/**
 * Copied and adapted from `struct mgaproblem` in `src/mga.h`.
 */
struct multiple_gravity_assist : mant::problem<double, 8> {
  multiple_gravity_assist() noexcept;

  std::array<const celestial_body*, 8> sequence;  // fly-by sequence
  std::array<int, 8> rev_flag;  // vector of flags for clockwise legs
  asteroid asteroid;
  double Isp;
  double mass;
  double DVlaunch;
};
}
