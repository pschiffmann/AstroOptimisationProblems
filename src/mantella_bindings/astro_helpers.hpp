#pragma once

#include <array>
#include <utility>

namespace multiple_gravity_assist {

typedef struct lambert_solution {
  std::array<double, 3> departure_velocity;
  std::array<double, 3> arrival_velocity;
} lambert_solution;

/**
 * This function an adapted version of function `LambertI` in file
 * `src/Lambert.h`. For more details and documentation, please take a look at
 * the original file.
 *
 * Throws `std::invalid_argument` if `t` is negative or zero.
 */
lambert_solution lambert(std::array<double, 3> r1_in,
                         std::array<double, 3> r2_in, double t, const double mu,
                         const int lw);

/**
 * This function an adapted version of function `PowSwingByInv` in file
 * `src/PowSwingByInv.h`.
 *
 * Returns a tuple (DV, rp).
 */
std::pair<double, double> PowSwingByInv(const double Vin, const double Vout,
                                        const double alpha);
}
