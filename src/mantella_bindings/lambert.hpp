#pragma once

#include <array>

namespace multiple_gravity_assist {

/**
 * This function an adapted version of function `LambertI` in file
 * `src/Lambert.h`. For more details and documentation, please take a look at
 * the original file.
 *
 * Throws `std::invalid_argument` if `t` is negative or zero.
 */
void lambert(std::array<double, 3> r1_in, std::array<double, 3> r2_in, double t,
             const double mu, const int lw,
             // OUTPUT
             std::array<double, 3>& v1, std::array<double, 3>& v2);
}
