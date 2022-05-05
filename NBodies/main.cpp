#include <stdio.h>

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>

#include "Integration/abmd/abmd.h"
#include "Utils/helper.h"
#include "Writer/writer.h"


void moonToBarycenter(std::vector<main_type> &rr, main_type *res)
{
  main_type earth_x = rr[earthNum * 6];
  main_type earth_y = rr[earthNum * 6 + 1];
  main_type earth_z = rr[earthNum * 6 + 2];
  main_type earth_vx = rr[earthNum * 6 + 3];
  main_type earth_vy = rr[earthNum * 6 + 4];
  main_type earth_vz = rr[earthNum * 6 + 5];

  main_type gcmoon_x = rr[moonNum * 6];
  main_type gcmoon_y = rr[moonNum * 6 + 1];
  main_type gcmoon_z = rr[moonNum * 6 + 2];
  main_type gcmoon_vx = rr[moonNum * 6 + 3];
  main_type gcmoon_vy = rr[moonNum * 6 + 4];
  main_type gcmoon_vz = rr[moonNum * 6 + 5];

  rr[moonNum * 6] = earth_x + gcmoon_x;
  rr[moonNum * 6 + 1] = earth_y + gcmoon_y;
  rr[moonNum * 6 + 2] = earth_z + gcmoon_z;
  rr[moonNum * 6 + 3] = earth_vx + gcmoon_vx;
  rr[moonNum * 6 + 4] = earth_vy + gcmoon_vy;
  rr[moonNum * 6 + 5] = earth_vz + gcmoon_vz;

  res[0] = gcmoon_x;
  res[1] = gcmoon_y;
  res[2] = gcmoon_z;
  res[3] = gcmoon_vx;
  res[4] = gcmoon_vy;
  res[5] = gcmoon_vz;
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

  // adjust initial positions and velocities

#ifdef SAVE_INV

  main_type *gcmoon = (main_type *)malloc(6 * sizeof(main_type));
  moonToBarycenter(rr, gcmoon);

  std::vector<main_type> pos_adj(3, 0);
  std::vector<main_type> vel_adj(3, 0);
  main_type sum_mu_star = 0;
  for(int i = 0; i < objects_counter; i++){
    main_type mu_star = 0;
    main_type mu_dot = 0;
    main_type sum = 0;
    main_type sum_dot = 0;

    for(int j = 0; j < objects_counter; j++){
      if(j != i){
        main_type dist = dist(rr[i*6], rr[i * 6 + 1], rr[i * 6 + 2],
                                      rr[j * 6], rr[j * 6 + 1], rr[j * 6 + 2]);
        sum_dot += masses[j] *  dotProduct(rr[i*6]   - rr[j*6], 
                                           rr[i*6+1] - rr[j*6+1], 
                                           rr[i*6+2] - rr[j*6+2],
                                           rr[i*6+3] + rr[j*6+3],
                                           rr[i*6+4] + rr[j*6+4],
                                           rr[i*6+5] + rr[j*6+5]) / (dist * dist * dist);
        sum += masses[j] / dist;
      }
    }
    mu_star = masses[i] * (1 + 1 / (2 * 299792.458 * 299792.458) * (vecLen2(rr[i * 6 + 3], rr[i * 6 + 4], rr[i * 6 + 5]) - sum));
    mu_dot = masses[i] * 1 / (2 * 299792.458 * 299792.458) * sum_dot;

    sum_mu_star += mu_star;
    pos_adj[0] += mu_star * rr[i * 6];
    pos_adj[1] += mu_star * rr[i * 6 + 1];
    pos_adj[2] += mu_star * rr[i * 6 + 2];

    vel_adj[0] += (mu_star * rr[i * 6 + 3] + mu_dot * rr[i * 6]);
    vel_adj[1] += (mu_star * rr[i * 6 + 4] + mu_dot * rr[i * 6 + 1]);
    vel_adj[2] += (mu_star * rr[i * 6 + 5] + mu_dot * rr[i * 6 + 2]);
  }

  for(int i = 0; i < objects_counter; i++){
    rr[i * 6]     -= pos_adj[0] / sum_mu_star;
    rr[i * 6 + 1] -= pos_adj[1] / sum_mu_star;
    rr[i * 6 + 2] -= pos_adj[2] / sum_mu_star;

    rr[i * 6 + 3] -= vel_adj[0] / sum_mu_star;
    rr[i * 6 + 4] -= vel_adj[1] / sum_mu_star;
    rr[i * 6 + 5] -= vel_adj[2] / sum_mu_star;
  }
  std::cout << pos_adj[0] << " " << pos_adj[1] << " " << pos_adj[2] << std::endl;
  std::cout << vel_adj[0] << " " << vel_adj[1] << " " << vel_adj[2] << std::endl;

  rr[moonNum * 6]     = gcmoon[0];
  rr[moonNum * 6 + 1] = gcmoon[1];
  rr[moonNum * 6 + 2] = gcmoon[2];
  rr[moonNum * 6 + 3] = gcmoon[3];
  rr[moonNum * 6 + 4] = gcmoon[4];
  rr[moonNum * 6 + 5] = gcmoon[5];

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
