#include "lambert.hpp"

#include <stdexcept>
#include "../Astro_Functions.h"
#include "vector3d_helpers.hpp"

void multiple_gravity_assist::lambert(std::array<double, 3> r1,
                                      std::array<double, 3> r2, double t,
                                      const double mu, const int lw,
                                      // OUTPUT
                                      std::array<double, 3>& v1,
                                      std::array<double, 3>& v2) {
  if (t <= 0) {
    throw std::invalid_argument(
        "ERROR in Lambert Solver: Negative Time in input.");
  }

  int i;
  double r2_vers[3];
  double ih_dum[3], ih[3], dum[3];

  double R = sqrt(dot_product(r1, r1));
  double V = sqrt(mu / R);
  double T = R / V;

  // working with non-dimensional radii and time-of-flight
  t /= T;
  r1 = div(r1, R);
  r2 = div(r2, R);

  // Evaluation of the relevant geometry parameters in non dimensional units
  // R2 module
  double r2_mod = sqrt(dot_product(r2, r2));

  double theta = acos(dot_product(r1, r2) / r2_mod);

  if (lw) theta = 2 * acos(-1.0) - theta;

  // non-dimensional chord
  double c = sqrt(1 + r2_mod * (r2_mod - 2.0 * cos(theta)));

  // non dimesnional semi-perimeter
  double s = (1 + r2_mod + c) / 2.0;

  // minimum energy ellipse semi major axis
  double am = s / 2.0;

  // lambda parameter defined in Battin's Book
  double lambda = sqrt(r2_mod) * cos(theta / 2.0) / s;

  double x;
  {
    // We start finding the log(x+1) value of the solution conic:
    // NO MULTI REV --> (1 SOL)
    //	inn1=-.5233;    //first guess point
    //  inn2=.5233;     //second guess point
    double x1 = log(0.4767);
    double x2 = log(1.5233);
    double y1 = log(x2tof(-.5233, s, c, lw)) - log(t);
    double y2 = log(x2tof(.5233, s, c, lw)) - log(t);

    // Increasing the tolerance does not bring any advantage as the
    // precision is usually greater anyway (due to the rectification of the tof
    // graph) except near particular cases such as parabolas in which cases a
    // lower precision allow for usual convergence.
    const double tolerance = 1e-11;
    double error = 1;

    // Newton iterations
    double x_new, y_new;
    while ((error > tolerance) && (y1 != y2)) {
      x_new = (x1 * y2 - y1 * x2) / (y2 - y1);
      y_new = logf(x2tof(expf(x_new) - 1, s, c, lw)) - logf(t);
      x1 = x2;
      y1 = y2;
      x2 = x_new;
      y2 = y_new;
      error = fabs(x1 - x_new);
    }
    x = expf(x_new) - 1;
  }

  // The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
  // now need the conic. As for transfer angles near to pi the lagrange
  // coefficient technique goes singular (dg approaches a zero/zero that is
  // numerically bad) we here use a different technique for those cases. When
  // the transfer angle is exactly equal to pi, then the ih unit vector is not
  // determined. The remaining equations, though, are still valid.

  double a = am / (1 - x * x);

  double eta, eta_squared;
  // psi evaluation
  if (x < 1) {
    // ellipse
    double beta = 2 * asin(sqrt((s - c) / (2 * a)));
    if (lw) beta = -beta;
    double alfa = 2 * acos(x);
    double psi = (alfa - beta) / 2;
    eta_squared = 2 * a * pow(sin(psi), 2) / s;
    eta = sqrt(eta_squared);
  } else {
    // hyperbola
    double beta = 2 * asinh(sqrt((c - s) / (2 * a)));
    if (lw) beta = -beta;
    double alfa = 2 * acosh(x);
    double psi = (alfa - beta) / 2;
    eta_squared = -2 * a * pow(sinh(psi), 2) / s;
    eta = sqrt(eta_squared);
  }

  // parameter of the solution
  double vr1, vt1, vt2, vr2;
  double p = (r2_mod / (am * eta_squared)) * pow(sin(theta / 2), 2);
  vr1 = (1 / (eta * sqrt(am))) * (2 * lambda * am - (lambda + x * eta));
  vett(r1.data(), r2.data(), ih_dum);
  vers(ih_dum, ih);

  if (lw) {
    for (i = 0; i < 3; i++) ih[i] = -ih[i];
  }

  vt1 = sqrt(p);
  vett(ih, r1.data(), dum);

  for (i = 0; i < 3; i++) v1[i] = vr1 * r1[i] + vt1 * dum[i];

  vt2 = vt1 / r2_mod;
  vr2 = -vr1 + (vt1 - vt2) / tan(theta / 2);
  vers(r2.data(), r2_vers);
  vett(ih, r2_vers, dum);
  for (i = 0; i < 3; i++) v2[i] = vr2 * r2[i] / r2_mod + vt2 * dum[i];

  v1 = mul(v1, V);
  v2 = mul(v2, V);
}
