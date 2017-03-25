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
void lambert(const double *r1_in, const double *r2_in, double t,
             const double mu,          // INPUT
             const int lw,             // INPUT
             double *v1, double *v2);  // OUTPUT
}
