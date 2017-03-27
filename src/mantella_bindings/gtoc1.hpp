#pragma once

#include <array>
#include <utility>

#include <mantella0>

namespace multiple_gravity_assist {
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

  double mu;
  double penalty;
  double penalty_coefficient;

  /**
   * This function is an adapted version of function
   * `Planet_Ephemerides_Analytical` in file `src/Pl_Eph_An.h`.
   *
   * Returns a tuple (position, velocity).
   */
  std::pair<std::array<double, 3>, std::array<double, 3>> ephemeris(
      const double mjd2000) const;

 private:
  celestial_body(double mu, double penalty,
                 double penalty_coefficient) noexcept;
};

struct asteroid {
 public:
  asteroid(std::array<double, 6> keplerian, double epoch, double mu) noexcept;

  std::array<double, 6> keplerian;
  double epoch;
  double mu;

  /**
   * This function is an adapted version of function `Custom_Eph` in file
   * `src/Pl_Eph_An.h`.
   *
   * Returns a tuple (position, velocity).
   */
  std::pair<std::array<double, 3>, std::array<double, 3>> ephemeris(
      const double jd) const;
};

/**
 * Copied and adapted from `struct mgaproblem` in `src/mga.h`.
 */
struct gtoc1 : mant::problem<double, 8> {
  gtoc1() noexcept;

  /**
   * fly-by sequence of planets. This sequence is 1 element smaller than the
   * problem dimension because the last section always targets the `asteroid`.
   */
  std::array<const celestial_body*, 7> sequence;
  std::array<bool, 8> rev_flag;  // vector of flags for clockwise legs
  asteroid asteroid;
  double Isp;
  double mass;
  double DVlaunch;
};
}
