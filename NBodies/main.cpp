#include <stdio.h>

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "Integration/abmd/abmd.h"
#include "Utils/helper.h"
#include "Writer/writer.h"

using namespace std;

int main() {
  std::ifstream infile("points.txt");
  long double x, y, z, vx, vy, vz, m;
  std::vector<main_type> rr;
  int objects_counter = 150;
  rr.resize(6 * objects_counter);
  std::vector<double> masses;
  masses.resize(objects_counter);
  std::string name;
  int k = 0;
  while ((infile >> name >> m >> x >> y >> z >> vx >> vy >> vz) && k < objects_counter) {
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
  //std::vector<main_type> moon_res(6, 0);
  //moonToGeocentric(rr.data(), moon_res.data());

  main_type total_mass = 0;
  
  // for (int i = 0; i < masses.size(); i++) {
  //   total_mass += masses[i];
  // }

  std::vector<main_type> init_center_mass(3, 0);

  // for (int i = 0; i < masses.size(); i++) {
  //   init_center_mass[0] += rr[i * 6] * masses[i];
  //   init_center_mass[1] += rr[i * 6 + 1] * masses[i];
  //   init_center_mass[2] += rr[i * 6 + 2] * masses[i];
  // }

  // init_center_mass[0] /= total_mass;
  // init_center_mass[1] /= total_mass;
  // init_center_mass[2] /= total_mass;

  std::vector<main_type> init_vel_mass(3, 0);

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

  main_type init_energy = 0;
  main_type init_k_energy = 0;
  // for (int i = 0; i < masses.size(); i++) {
  //   init_k_energy +=
  //       masses[i] *
  //       (rr[i * 6 + 3] * rr[i * 6 + 3] + rr[i * 6 + 4] * rr[i * 6 + 4] +
  //        rr[i * 6 + 5] * rr[i * 6 + 5]) /
  //       2.0;
  // }
  main_type init_p_energy = 0;
  // for (int i = 0; i < masses.size(); i++) {
  //   for (int j = i + 1; j < masses.size(); j++) {
  //     init_p_energy += masses[i] * masses[j] /
  //                      sqrt((rr[i * 6] - rr[j * 6]) * (rr[i * 6] - rr[j * 6])
  //                      +
  //                           (rr[i * 6 + 1] - rr[j * 6 + 1]) *
  //                               (rr[i * 6 + 1] - rr[j * 6 + 1]) +
  //                           (rr[i * 6 + 2] - rr[j * 6 + 2]) *
  //                               (rr[i * 6 + 2] - rr[j * 6 + 2]));
  //   }
  // }
  // init_energy = init_k_energy - init_p_energy;

  std::vector<main_type> init_impulse_moment(3, 0);
  // for (int i = 0; i < masses.size(); i++) {
  //   main_type temp_x = rr[i * 6 + 3] * masses[i];
  //   main_type temp_y = rr[i * 6 + 4] * masses[i];
  //   main_type temp_z = rr[i * 6 + 5] * masses[i];

  //   init_impulse_moment[0] += rr[i * 6 + 1] * temp_z - rr[i * 6 + 2] *
  //   temp_y; init_impulse_moment[1] += rr[i * 6 + 2] * temp_x - rr[i * 6] *
  //   temp_z; init_impulse_moment[2] += rr[i * 6] * temp_y - rr[i * 6 + 1] *
  //   temp_x;
  // }

  // rr[moonNum * 6] = moon_res[0];
  // rr[moonNum * 6 + 1] = moon_res[1];
  // rr[moonNum * 6 + 2] = moon_res[2];
  // rr[moonNum * 6 + 3] = moon_res[3];
  // rr[moonNum * 6 + 4] = moon_res[4];
  // rr[moonNum * 6 + 5] = moon_res[5];

  double h = 0.03125;
  auto start = std::chrono::high_resolution_clock::now();
  double t = 365 * YEARS;
  int sol_size = 2 * (int)(1 + t / h) - 1;

  main_type *diff = (main_type *)malloc(sizeof(main_type) * sol_size);

  ABMD_calc_diff(rr, masses, h, t, init_energy, init_impulse_moment.data(),
                 init_center_mass.data(), diff);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "time: " << duration.count() << " microseconds \n";

  return 0;
}
