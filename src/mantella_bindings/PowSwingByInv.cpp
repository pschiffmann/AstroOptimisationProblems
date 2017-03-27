#include "PowSwingByInv.hpp"

#include <math.h>

std::pair<double, double> multiple_gravity_assist::PowSwingByInv(
    const double Vin, const double Vout, const double alpha) {
  double DV, rp;
  const int maxiter = 30;
  int i = 0;
  double err = 1.0;
  double f, df;  // function and its derivative
  double rp_new;
  const double tolerance = 1e-8;

  double aIN = 1.0 / pow(Vin, 2);    // semimajor axis of the incoming hyperbola
  double aOUT = 1.0 / pow(Vout, 2);  // semimajor axis of the incoming hyperbola

  rp = 1.0;
  while ((err > tolerance) && (i < maxiter)) {
    i++;
    f = asin(aIN / (aIN + rp)) + asin(aOUT / (aOUT + rp)) - alpha;
    df = -aIN / sqrt((rp + 2 * aIN) * rp) / (aIN + rp) -
         aOUT / sqrt((rp + 2 * aOUT) * rp) / (aOUT + rp);
    rp_new = rp - f / df;
    if (rp_new > 0) {
      err = fabs(rp_new - rp);
      rp = rp_new;
    } else
      rp /= 2.0;
  }

  // Evaluation of the DV
  DV = fabs(sqrt(Vout * Vout + (2.0 / rp)) - sqrt(Vin * Vin + (2.0 / rp)));

  return {DV, rp};
}
