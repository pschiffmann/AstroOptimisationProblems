#include "gtoc1.hpp"

#include <iterator>
#include <numeric>

#include "../Pl_Eph_An.h"
#include "../PowSwingByInv.h"
#include "lambert.hpp"
#include "vector3d_helpers.hpp"

using namespace multiple_gravity_assist;

const celestial_body celestial_body::SUN(1.32712428e11, 0, 0, 0);
const celestial_body celestial_body::MERCURY(22321, 0, 0, 1);
const celestial_body celestial_body::VENUS(324860, 6351.8, 0.01, 2);
const celestial_body celestial_body::EARTH(398601.19, 6778.1, 0.01, 3);
const celestial_body celestial_body::MARS(42828.3, 6000, 0.01, 4);
const celestial_body celestial_body::JUPITER(126.7e6, 600000, 0.001, 5);
const celestial_body celestial_body::SATURN(37.9e6, 70000, 0.01, 6);
const celestial_body celestial_body::URANUS(5.78e6, 0, 0, 7);
const celestial_body celestial_body::NEPTUNE(6.8e6, 0, 0, 8);

celestial_body::celestial_body(double mu, double penalty,
                               double penalty_coefficient,
                               unsigned char distance_from_sun) noexcept
    : mu(mu),
      penalty(penalty),
      penalty_coefficient(penalty_coefficient),
      distance_from_sun(distance_from_sun) {}

asteroid::asteroid(std::array<double, 6> keplerian, double epoch,
                   double mu) noexcept
    : keplerian(keplerian), epoch(epoch), mu(mu) {}

gtoc1::gtoc1() noexcept
    : sequence({&celestial_body::EARTH, &celestial_body::VENUS,
                &celestial_body::EARTH, &celestial_body::VENUS,
                &celestial_body::EARTH, &celestial_body::JUPITER,
                &celestial_body::SATURN}),
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
  /**
   * Adapted version of function `MGA` from `src/mga.cpp`.
   */
  this->objective_function = [this](const array<double, 8> &parameter) {
    std::array<double, 6> rp{};
    std::array<double, 8> DV{};
    const int n = 8;

    const double g = 9.80665 / 1000.0;  // Gravity

    // r and  v in heliocentric coordinate system
    // position
    std::array<std::array<double, 3>, 8> r;
    // velocity
    std::array<std::array<double, 3>, 8> v;

    {
      double totalTime = 0;
      for (size_t i = 0; i < 7; i++) {
        totalTime += parameter[i];
        Planet_Ephemerides_Analytical(totalTime,
                                      sequence.at(i)->distance_from_sun,
                                      r[i].data(), v[i].data());
      }
      totalTime += parameter[7];
      Custom_Eph(totalTime + 2451544.5, asteroid.epoch,
                 asteroid.keplerian.data(), r[7].data(), v[7].data());
    }

    std::array<double, 3> current_section_departure_velocity;
    std::array<double, 3> current_section_arrival_velocity;
    for (size_t i = 0; i <= n - 2; i++) {
      std::array<double, 3> last_section_departure_velocity =
          current_section_departure_velocity;
      std::array<double, 3> last_section_arrival_velocity =
          current_section_arrival_velocity;

      bool longWay =
          cross_product(r[i], r[i + 1])[2] > 0 ? rev_flag[i] : !rev_flag[i];

      lambert(r[i], r[i + 1], parameter[i + 1] * 24 * 60 * 60,
              celestial_body::SUN.mu, longWay,
              // OUTPUT
              current_section_departure_velocity.data(),
              current_section_arrival_velocity.data());

      if (i == 0) {
        // Earth launch
        DV[0] = norm(sub(current_section_departure_velocity, v[0]));
      } else {
        double Vin = norm(sub(last_section_arrival_velocity, v[i]));
        double Vout = norm(sub(current_section_departure_velocity, v[i]));

        // calculation of delta V at pericenter
        PowSwingByInv(
            Vin, Vout,
            acos(dot_product(sub(last_section_arrival_velocity, v[i]),
                             sub(current_section_departure_velocity, v[i])) /
                 (Vin * Vout)),
            DV[i], rp[i - 1]);
        rp[i - 1] *= sequence[i]->mu;
      }
    }

    double DVtot = std::accumulate<double *, double>(std::next(DV.begin()),
                                                     std::prev(DV.end()), 0);

    // Build Penalty
    for (size_t i = 0; i < n - 2; i++) {
      if (rp[i] < sequence[i + 1]->penalty)
        DVtot += sequence[i + 1]->penalty_coefficient *
                 fabs(rp[i] - sequence[i + 1]->penalty);
    }

    // Launcher Constraint
    if (DV[0] > DVlaunch) DVtot += (DV[0] - DVlaunch);

    // Evaluation of satellite final mass
    double final_mass = mass * exp(-DVtot / (Isp * g));

    // arrival relative velocity at the asteroid;
    auto arrival_velocity = sub(v[n - 1], current_section_arrival_velocity);

    return -final_mass * fabs(dot_product(arrival_velocity, v[n - 1]));
  };
}
