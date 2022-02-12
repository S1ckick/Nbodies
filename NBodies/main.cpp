#include <stdio.h>

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

//#define NUMBER_DOUBLE 1
#define NUMBER_DOUBLE_DOUBLE 1

#ifdef NUMBER_DOUBLE_DOUBLE
#include <qd/dd_real.h>
#include <qd/fpu.h>

using current_type = dd_real;
#else
using current_type = long double;
#endif

#include "Integration/abmd/abmd.h"
#include "Utils/helper.h"
#include "Writer/writer.h"

using namespace std;
// using current_type = double;

int main() {
  std::ifstream infile("../NBodies/points.txt");
  long double x, y, z, vx, vy, vz, m;
  std::vector<current_type> rr;
  rr.resize(6 * 16);
  std::vector<double> masses;
  masses.resize(16);
  std::string name;
  int k = 0;
  while (infile >> name >> m >> x >> y >> z >> vx >> vy >> vz) {
    masses[k] = m;

    rr[k * 6] = x;
    rr[k * 6 + 1] = y;
    rr[k * 6 + 2] = z;

    rr[k * 6 + 3] = vx;
    rr[k * 6 + 4] = vy;
    rr[k * 6 + 5] = vz;
    std::cout << "parsed " << name << " : " << x << " " << y << " " << z << " "
              << vx << " " << vy << " " << vz << " " << m << std::endl;
    k++;
  }

  /// for energy fix need change moon center
  //std::vector<current_type> moon_res(6, 0);
  //moonToGeocentric(rr.data(), moon_res.data());

  current_type total_mass = 0;
  
  // for (int i = 0; i < masses.size(); i++) {
  //   total_mass += masses[i];
  // }

  std::vector<current_type> init_center_mass(3, 0);

  // for (int i = 0; i < masses.size(); i++) {
  //   init_center_mass[0] += rr[i * 6] * masses[i];
  //   init_center_mass[1] += rr[i * 6 + 1] * masses[i];
  //   init_center_mass[2] += rr[i * 6 + 2] * masses[i];
  // }

  // init_center_mass[0] /= total_mass;
  // init_center_mass[1] /= total_mass;
  // init_center_mass[2] /= total_mass;

  std::vector<current_type> init_vel_mass(3, 0);

  // for (int i = 0; i < masses.size(); i++) {
  //   init_vel_mass[0] += rr[i * 6 + 3] * masses[i];
  //   init_vel_mass[1] += rr[i * 6 + 4] * masses[i];
  //   init_vel_mass[2] += rr[i * 6 + 5] * masses[i];
  // }

  // init_vel_mass[0] /= total_mass;
  // init_vel_mass[1] /= total_mass;
  // init_vel_mass[2] /= total_mass;

  // for (int i = 0; i < masses.size(); i++) {
  //   rr[i * 6] -= init_center_mass[0];
  //   rr[i * 6 + 1] -= init_center_mass[1];
  //   rr[i * 6 + 2] -= init_center_mass[2];
  //   rr[i * 6 + 3] -= init_vel_mass[0];
  //   rr[i * 6 + 4] -= init_vel_mass[1];
  //   rr[i * 6 + 5] -= init_vel_mass[2];
  // }

  current_type init_energy = 0;
  current_type init_k_energy = 0;
  // for (int i = 0; i < masses.size(); i++) {
  //   init_k_energy +=
  //       masses[i] *
  //       (rr[i * 6 + 3] * rr[i * 6 + 3] + rr[i * 6 + 4] * rr[i * 6 + 4] +
  //        rr[i * 6 + 5] * rr[i * 6 + 5]) /
  //       2.0;
  // }
  current_type init_p_energy = 0;
  // for (int i = 0; i < masses.size(); i++) {
  //   for (int j = i + 1; j < masses.size(); j++) {
  //     init_p_energy += masses[i] * masses[j] /
  //                      sqrt((rr[i * 6] - rr[j * 6]) * (rr[i * 6] - rr[j * 6]) +
  //                           (rr[i * 6 + 1] - rr[j * 6 + 1]) *
  //                               (rr[i * 6 + 1] - rr[j * 6 + 1]) +
  //                           (rr[i * 6 + 2] - rr[j * 6 + 2]) *
  //                               (rr[i * 6 + 2] - rr[j * 6 + 2]));
  //   }
  // }
  // init_energy = init_k_energy - init_p_energy;

  std::vector<current_type> init_impulse_moment(3, 0);
  // for (int i = 0; i < masses.size(); i++) {
  //   current_type temp_x = rr[i * 6 + 3] * masses[i];
  //   current_type temp_y = rr[i * 6 + 4] * masses[i];
  //   current_type temp_z = rr[i * 6 + 5] * masses[i];

  //   init_impulse_moment[0] += rr[i * 6 + 1] * temp_z - rr[i * 6 + 2] * temp_y;
  //   init_impulse_moment[1] += rr[i * 6 + 2] * temp_x - rr[i * 6] * temp_z;
  //   init_impulse_moment[2] += rr[i * 6] * temp_y - rr[i * 6 + 1] * temp_x;
  // }

  // rr[moonNum * 6] = moon_res[0];
  // rr[moonNum * 6 + 1] = moon_res[1];
  // rr[moonNum * 6 + 2] = moon_res[2];
  // rr[moonNum * 6 + 3] = moon_res[3];
  // rr[moonNum * 6 + 4] = moon_res[4];
  // rr[moonNum * 6 + 5] = moon_res[5];

  double h = 0.03125;
  auto start = std::chrono::high_resolution_clock::now();
  double t = 365 * 40;
  int sol_size = 2 * (int)(1 + t / h) - 1;

  current_type *diff = (current_type *)malloc(sizeof(current_type) * sol_size);

  ABMD_calc_diff(rr, masses, h, t, init_energy, init_impulse_moment.data(),
                 init_center_mass.data(), diff);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "time: " << duration.count() << " microseconds \n";

  return 0;
}
