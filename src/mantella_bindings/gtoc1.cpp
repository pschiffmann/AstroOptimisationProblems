#include "gtoc1.hpp"

using namespace gtoc1;

const celestial_body celestial_body::SUN(1.32712428e11, 0, 0);
const celestial_body celestial_body::MERCURY(22321, 0, 0);
const celestial_body celestial_body::VENUS(324860, 6351.8, 0.01);
const celestial_body celestial_body::EARTH(398601.19, 6778.1, 0.01);
const celestial_body celestial_body::MARS(42828.3, 6000, 0.01);
const celestial_body celestial_body::JUPITER(126.7e6, 600000, 0.001);
const celestial_body celestial_body::SATURN(37.9e6, 70000, 0.01);
const celestial_body celestial_body::URANUS(5.78e6, 0, 0);
const celestial_body celestial_body::NEPTUNE(6.8e6, 0, 0);

celestial_body::celestial_body(double mu, double penalty,
                               double penalty_coefficient) noexcept
    : mu(mu), penalty(penalty), penalty_coefficient(penalty_coefficient) {}

asteroid::asteroid(std::array<double, 6> keplerian, double epoch,
                   double mu) noexcept
    : keplerian(keplerian), epoch(epoch), mu(mu) {}

multiple_gravity_assist::multiple_gravity_assist() noexcept
    : sequence({&celestial_body::EARTH, &celestial_body::VENUS,
                &celestial_body::EARTH, &celestial_body::VENUS,
                &celestial_body::EARTH, &celestial_body::JUPITER,
                &celestial_body::SATURN,
                &celestial_body::NEPTUNE /*REPLACE WITH asteroid*/}),
      rev_flag({0, 0, 0, 0, 0, 0, 1, 0}),
      asteroid(
          {{2.5897261, 0.2734625, 6.40734, 128.34711, 264.78691, 320.479555},
           53600.0,
           0.0}),
      Isp(2500.0),
      mass(1500.0),
      DVlaunch(2.5),
      problem() {
  // this->lower_bounds = TODO
  // this->upper_bounds = TODO
  this->objective_function = [this](const auto &parameter) {
    std::array<double, 6> rp;
    std::vector<double> Delta_V(8);
    double obj = 0;
    // MGA(x, problem, rp, Delta_V, obj);
    return obj;
  };
}
