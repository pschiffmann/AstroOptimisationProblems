#pragma once

#include <utility>

namespace multiple_gravity_assist {
/**
 * This function an adapted version of function `PowSwingByInv` in file
 * `src/PowSwingByInv.h`.
 *
 * Returns a tuple (DV, rp).
 */
std::pair<double, double> PowSwingByInv(const double Vin, const double Vout,
                                        const double alpha);
}
