#include "gtoc1.hpp"

#include <math.h>
#include <numeric>

#include "../Lambert.h"
#include "../Pl_Eph_An.h"
#include "../PowSwingByInv.h"

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
    std::array<double, 6> rp;
    std::array<double, 8> DV;
    double obj_funct;
    const int n = 8;

    double DVtot = 0;
    double Dum_Vec[3], Vin, Vout;
    double V_Lamb[2][2][3], dot_prod;
    double a, p, theta, alfa;
    double DVrel, DVarr = 0;

    // only used for asteroid impact (ex: gtoc1)
    const double initial_mass = mass;   // Satellite initial mass [Kg]
    double final_mass;                  // satelite final mass
    const double g = 9.80665 / 1000.0;  // Gravity

    double *vec, *rec;
    vector<double *> r;  // {0...n-1} position
    vector<double *> v;  // {0...n-1} velocity

    double T = 0.0;  // total time

    int i_count, j_count, lw;

    int iter = 0;

    {
      for (i_count = 0; i_count < n; i_count++) {
        vec = new double[3];  // velocity and position are 3 D vector
        rec = new double[3];
        r.push_back(vec);
        v.push_back(rec);

        DV[i_count] = 0.0;
      }

      T = 0;
      for (i_count = 0; i_count < 7; i_count++) {
        T += parameter[i_count];
        Planet_Ephemerides_Analytical(
            T, sequence.at(i_count)->distance_from_sun, r[i_count],
            v[i_count]);  // r and  v in heliocentric coordinate system
      }
      T += parameter[7];
      Custom_Eph(T + 2451544.5, asteroid.epoch, asteroid.keplerian.data(), r[7],
                 v[7]);

      vett(r[0], r[1], Dum_Vec);

      if (Dum_Vec[2] > 0)
        lw = (rev_flag[0] == 0) ? 0 : 1;
      else
        lw = (rev_flag[0] == 0) ? 1 : 0;

      LambertI(r[0], r[1], parameter[1] * 24 * 60 * 60, celestial_body::SUN.mu,
               lw,                                              // INPUT
               V_Lamb[0][0], V_Lamb[0][1], a, p, theta, iter);  // OUTPUT
      DV[0] = norm(V_Lamb[0][0], v[0]);                         // Earth launch

      for (i_count = 1; i_count <= n - 2; i_count++) {
        vett(r[i_count], r[i_count + 1], Dum_Vec);

        if (Dum_Vec[2] > 0)
          lw = (rev_flag[i_count] == 0) ? 0 : 1;
        else
          lw = (rev_flag[i_count] == 0) ? 1 : 0;

        /*if (i_count%2 != 0)	{*/
        LambertI(r[i_count], r[i_count + 1],
                 parameter[i_count + 1] * 24 * 60 * 60, celestial_body::SUN.mu,
                 lw,                                              // INPUT
                 V_Lamb[1][0], V_Lamb[1][1], a, p, theta, iter);  // OUTPUT

        // norm first perform the subtraction of vet1-vet2 and the evaluate
        // ||...||
        Vin = norm(V_Lamb[0][1], v[i_count]);
        Vout = norm(V_Lamb[1][0], v[i_count]);

        dot_prod = 0.0;
        for (int i = 0; i < 3; i++) {
          dot_prod += (V_Lamb[0][1][i] - v[i_count][i]) *
                      (V_Lamb[1][0][i] - v[i_count][i]);
        }
        alfa = acos(dot_prod / (Vin * Vout));

        // calculation of delta V at pericenter
        PowSwingByInv(Vin, Vout, alfa, DV[i_count], rp[i_count - 1]);

        rp[i_count - 1] *= sequence[i_count]->mu;

        if (i_count != n - 2)  // swap
          for (j_count = 0; j_count < 3; j_count++) {
            V_Lamb[0][0][j_count] = V_Lamb[1][0][j_count];  // [j_count];
            V_Lamb[0][1][j_count] = V_Lamb[1][1][j_count];  // [j_count];
          }
      }
    }

    for (i_count = 0; i_count < 3; i_count++)
      Dum_Vec[i_count] = v[n - 1][i_count] - V_Lamb[1][1][i_count];

    DVrel = norm2(Dum_Vec);

    DVarr = DVrel;

    DVtot = 0;

    for (i_count = 1; i_count < n - 1; i_count++) DVtot += DV[i_count];

    // Build Penalty
    for (i_count = 0; i_count < n - 2; i_count++)
      if (rp[i_count] < sequence[i_count + 1]->penalty)
        DVtot += sequence[i_count + 1]->penalty_coefficient *
                 fabs(rp[i_count] - sequence[i_count + 1]->penalty);

    // Launcher Constraint
    if (DV[0] > DVlaunch) DVtot += (DV[0] - DVlaunch);

    // Evaluation of satellite final mass
    obj_funct = final_mass = initial_mass * exp(-DVtot / (Isp * g));

    // V asteroid - V satellite
    for (i_count = 0; i_count < 3; i_count++)
      Dum_Vec[i_count] =
          v[n - 1][i_count] -
          V_Lamb[1][1][i_count];  // arrival relative velocity at the asteroid;

    dot_prod = 0;
    for (i_count = 0; i_count < 3; i_count++)
      dot_prod += Dum_Vec[i_count] * v[n - 1][i_count];

    obj_funct = -(final_mass)*fabs(dot_prod);

    // final clean
    for (i_count = 0; i_count < n; i_count++) {
      delete[] r[i_count];
      delete[] v[i_count];
    }
    r.clear();
    v.clear();
    return obj_funct;
  };
}