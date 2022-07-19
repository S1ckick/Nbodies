#include <stdio.h>

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "Integration/abmd/abmd.h"
#include "Utils/helper.h"
#include "Writer/writer.h"

void moonToBarycenter(std::vector<main_type> &rr, std::vector<main_type> &res) {
  res[0] = rr[moonNum * 6];
  res[1] = rr[moonNum * 6 + 1];
  res[2] = rr[moonNum * 6 + 2];
  res[3] = rr[moonNum * 6 + 3];
  res[4] = rr[moonNum * 6 + 4];
  res[5] = rr[moonNum * 6 + 5];

  rr[moonNum * 6] += rr[earthNum * 6];
  rr[moonNum * 6 + 1] += rr[earthNum * 6 + 1];
  rr[moonNum * 6 + 2] += rr[earthNum * 6 + 2];
  rr[moonNum * 6 + 3] += rr[earthNum * 6 + 3];
  rr[moonNum * 6 + 4] += rr[earthNum * 6 + 4];
  rr[moonNum * 6 + 5] += rr[earthNum * 6 + 5];
}

void transform_x_to_center(std::vector<double> &masses,
                           std::vector<main_type> &rr, int n) {
  std::vector<main_type> gcmoon(6, 0);
  moonToBarycenter(rr, gcmoon);
#ifdef RELATIVISTIC
#ifdef NUMBER_DOUBLE_DOUBLE
  main_type revLS2 = main_type("3.335661199676477669820576657401368279e-05");
#else
#if NUMBER_DOUBLE == 1
  main_type revLS2 = 3.335661199676477669820577e-05L;
#else
  main_type revLS2 = 3.335661199676477670e-005;
#endif
#endif
  for (int ii = 0; ii < n; ii++) {
    std::vector<main_type> pos_adj(3, 0);
    std::vector<main_type> vel_adj(3, 0);
    main_type sum_mu_star = 0;

    for (int i = 0; i < masses.size(); i++) {
      main_type mu_star = 0;
      main_type mu_dot = 0;
      main_type sum = 0;
      main_type sum_dot = 0;

      int jn = i < barrier ? masses.size() : barrier;
      for (int j = 0; j < jn; j++) {
        if (j != i) {
          main_type dist = dist(rr[j * 6], rr[j * 6 + 1], rr[j * 6 + 2],
                                rr[i * 6], rr[i * 6 + 1], rr[i * 6 + 2]);
          sum_dot +=
              ((main_type)masses[j]) *
              dotProduct(rr[j * 6] - rr[i * 6], rr[j * 6 + 1] - rr[i * 6 + 1],
                         rr[j * 6 + 2] - rr[i * 6 + 2],
                         rr[j * 6 + 3] + rr[i * 6 + 3],
                         rr[j * 6 + 4] + rr[i * 6 + 4],
                         rr[j * 6 + 5] + rr[i * 6 + 5]) /
              (dist * dist * dist);
          sum += ((main_type)masses[j]) / dist;
        }
      }
      mu_star =
          ((main_type)masses[i]) *
          (1. +
           0.5 * revLS2 *
               (vecLen2(rr[i * 6 + 3], rr[i * 6 + 4], rr[i * 6 + 5]) - sum));
      mu_dot = ((main_type)masses[i]) * 0.5 * revLS2 * sum_dot;

      sum_mu_star += mu_star;
      pos_adj[0] += mu_star * rr[i * 6];
      pos_adj[1] += mu_star * rr[i * 6 + 1];
      pos_adj[2] += mu_star * rr[i * 6 + 2];

      vel_adj[0] += (mu_star * rr[i * 6 + 3] + mu_dot * rr[i * 6]);
      vel_adj[1] += (mu_star * rr[i * 6 + 4] + mu_dot * rr[i * 6 + 1]);
      vel_adj[2] += (mu_star * rr[i * 6 + 5] + mu_dot * rr[i * 6 + 2]);
    }

    // for (int i = 0; i < masses.size(); i++) {
    //   rr[i * 6] -= pos_adj[0] / sum_mu_star;
    //   rr[i * 6 + 1] -= pos_adj[1] / sum_mu_star;
    //   rr[i * 6 + 2] -= pos_adj[2] / sum_mu_star;

    //   rr[i * 6 + 3] -= vel_adj[0] / sum_mu_star;
    //   rr[i * 6 + 4] -= vel_adj[1] / sum_mu_star;
    //   rr[i * 6 + 5] -= vel_adj[2] / sum_mu_star;
    // }
    std::cout << "Barycentre pos: " << pos_adj[0] << " " << pos_adj[1] << " "
              << pos_adj[2] << std::endl;
    std::cout << "Barycentre vel" << vel_adj[0] << " " << vel_adj[1] << " "
              << vel_adj[2] << std::endl;
  }
#else
  std::vector<main_type> pos_adj(3, 0);
  std::vector<main_type> vel_adj(3, 0);
  main_type total_mass = 0;
  for (int i = 0; i < masses.size(); i++) {
    total_mass += masses[i];
    pos_adj[0] += rr[i * 6] * masses[i];
    pos_adj[1] += rr[i * 6 + 1] * masses[i];
    pos_adj[2] += rr[i * 6 + 2] * masses[i];

    vel_adj[0] += rr[i * 6 + 3] * masses[i];
    vel_adj[1] += rr[i * 6 + 4] * masses[i];
    vel_adj[2] += rr[i * 6 + 5] * masses[i];
  }

  pos_adj[0] /= total_mass;
  pos_adj[1] /= total_mass;
  pos_adj[2] /= total_mass;
  vel_adj[0] /= total_mass;
  vel_adj[1] /= total_mass;
  vel_adj[2] /= total_mass;

  for (int i = 0; i < masses.size(); i++) {
    rr[i * 6] -= pos_adj[0];
    rr[i * 6 + 1] -= pos_adj[1];
    rr[i * 6 + 2] -= pos_adj[2];

    rr[i * 6 + 3] -= vel_adj[0];
    rr[i * 6 + 4] -= vel_adj[1];
    rr[i * 6 + 5] -= vel_adj[2];
  }

#endif
  rr[moonNum * 6] = gcmoon[0];
  rr[moonNum * 6 + 1] = gcmoon[1];
  rr[moonNum * 6 + 2] = gcmoon[2];
  rr[moonNum * 6 + 3] = gcmoon[3];
  rr[moonNum * 6 + 4] = gcmoon[4];
  rr[moonNum * 6 + 5] = gcmoon[5];
}

using namespace std;

int main() {
  std::ifstream infile("pointmasses.txt");
  long double x, y, z, vx, vy, vz, m;
  std::vector<main_type> rr;
  int objects_counter = 660;
  rr.resize(6 * objects_counter);
  std::vector<double> masses;
  masses.resize(objects_counter);
  std::string name;
  int k = 0;
  while ((infile >> name >> m >> x >> y >> z >> vx >> vy >> vz) &&
         k < objects_counter) {
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

#ifdef SAVE_INV
  // adjust initial positions and velocities
  transform_x_to_center(masses, rr, 3);
  // END of adjustment
#endif

  double h = 0.03125;
  auto start = std::chrono::high_resolution_clock::now();
  double t = 365 * YEARS;

  ABMD_calc_diff(rr, masses, h, t);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "time: " << duration.count() << " microseconds \n";

  return 0;
}
