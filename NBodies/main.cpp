#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#define NUMBER_DOUBLE 1
//#define NUMBER_DOUBLE_DOUBLE 1

#include "Integration/methods.h"
#include "Writer/writer.h"

#ifdef NUMBER_DOUBLE_DOUBLE
#include <qd/dd_real.h>
#include <qd/fpu.h>

using current_type = dd_real;
#else
using current_type = long double;
#endif

#include "Utils/helper.h"

using namespace std;
using current_type = long double;

void moonToGeocentric(std::vector<current_type> &rr,
                      std::vector<current_type> &res) {
  current_type earth_x = rr[earthNum * 6];
  current_type earth_y = rr[earthNum * 6 + 1];
  current_type earth_z = rr[earthNum * 6 + 2];
  current_type earth_vx = rr[earthNum * 6 + 3];
  current_type earth_vy = rr[earthNum * 6 + 4];
  current_type earth_vz = rr[earthNum * 6 + 5];

  current_type gcmoon_x = rr[moonNum * 6];
  current_type gcmoon_y = rr[moonNum * 6 + 1];
  current_type gcmoon_z = rr[moonNum * 6 + 2];
  current_type gcmoon_vx = rr[moonNum * 6 + 3];
  current_type gcmoon_vy = rr[moonNum * 6 + 4];
  current_type gcmoon_vz = rr[moonNum * 6 + 5];

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

int main() {
#ifdef NUMBER_DOUBLE_DOUBLE
  unsigned int oldcw;
  fpu_fix_start(&oldcw);
#endif

  std::ifstream infile("../NBodies/points.txt");
  long double x, y, z, vx, vy, vz, m;
  std::vector<current_type> rr;
  rr.resize(6 * 16);
  std::vector<current_type> masses;
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
  std::vector<current_type> moon_res(6, 0);
  moonToGeocentric(rr, moon_res);

  current_type total_mass = 0;
  for (int i = 0; i < masses.size(); i++) {
    total_mass += masses[i];
  }

  std::vector<current_type> init_center_mass(3, 0);

  for (int i = 0; i < masses.size(); i++) {
    init_center_mass[0] += rr[i * 6] * masses[i];
    init_center_mass[1] += rr[i * 6 + 1] * masses[i];
    init_center_mass[2] += rr[i * 6 + 2] * masses[i];
  }

  init_center_mass[0] /= total_mass;
  init_center_mass[1] /= total_mass;
  init_center_mass[2] /= total_mass;

  std::vector<current_type> init_vel_mass(3, 0);

  for (int i = 0; i < masses.size(); i++) {
    init_vel_mass[0] += rr[i * 6 + 3] * masses[i];
    init_vel_mass[1] += rr[i * 6 + 4] * masses[i];
    init_vel_mass[2] += rr[i * 6 + 5] * masses[i];
  }

  init_vel_mass[0] /= total_mass;
  init_vel_mass[1] /= total_mass;
  init_vel_mass[2] /= total_mass;

  for (int i = 0; i < masses.size(); i++) {
    rr[i * 6] -= init_center_mass[0];
    rr[i * 6 + 1] -= init_center_mass[1];
    rr[i * 6 + 2] -= init_center_mass[2];
    rr[i * 6 + 3] -= init_vel_mass[3];
    rr[i * 6 + 4] -= init_vel_mass[4];
    rr[i * 6 + 5] -= init_vel_mass[5];
  }

  current_type init_energy = 0;
  current_type init_k_energy = 0;
  for (int i = 0; i < masses.size(); i++) {
    init_k_energy +=
        masses[i] *
        (rr[i * 6 + 3] * rr[i * 6 + 3] + rr[i * 6 + 4] * rr[i * 6 + 4] +
         rr[i * 6 + 5] * rr[i * 6 + 5]) /
        2.0;
  }
  current_type init_p_energy = 0;
  for (int i = 0; i < masses.size(); i++) {
    for (int j = i + 1; j < masses.size(); j++) {
      init_p_energy +=
          masses[i] * masses[j] /
          std::sqrt((rr[i * 6] - rr[j * 6]) * (rr[i * 6] - rr[j * 6]) +
                     (rr[i * 6 + 1] - rr[j * 6 + 1]) *
                         (rr[i * 6 + 1] - rr[j * 6 + 1]) +
                     (rr[i * 6 + 2] - rr[j * 6 + 2]) *
                         (rr[i * 6 + 2] - rr[j * 6 + 2]));
    }
  }
  init_energy = init_k_energy - init_p_energy;

  std::vector<current_type> init_impulse_moment(3, 0);
  for (int i = 0; i < masses.size(); i++) {
    current_type temp_x = rr[i * 6 + 3] * masses[i];
    current_type temp_y = rr[i * 6 + 4] * masses[i];
    current_type temp_z = rr[i * 6 + 5] * masses[i];

    init_impulse_moment[0] += rr[i * 6 + 1] * temp_z - rr[i * 6 + 2] * temp_y;
    init_impulse_moment[1] += rr[i * 6 + 2] * temp_x - rr[i * 6] * temp_z;
    init_impulse_moment[2] += rr[i * 6] * temp_y - rr[i * 6 + 1] * temp_x;
  }

  rr[moonNum * 6] = moon_res[0];
  rr[moonNum * 6 + 1] = moon_res[1];
  rr[moonNum * 6 + 2] = moon_res[2];
  rr[moonNum * 6 + 3] = moon_res[3];
  rr[moonNum * 6 + 4] = moon_res[4];
  rr[moonNum * 6 + 5] = moon_res[5];

  std::vector<current_type> data_bodies, data_energy, data_impulse_moment,
      data_center;

  current_type h(0.0625);
  int iterations = 5000;

  data_bodies.resize(masses.size() * iterations * 3);
  data_energy.resize(iterations);
  data_impulse_moment.resize(iterations);
  data_center.resize(iterations);

  std::vector<current_type> coefs = initDDCoef<current_type>();
  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < iterations; i++) {
    RungeKutta4(rr, masses, h);

    for (int j = 0; j < masses.size(); j++) {
      data_bodies[3 * (j + masses.size() * i)] = rr[j * 6];
      data_bodies[1 + 3 * (j + masses.size() * i)] = rr[j * 6 + 1];
      data_bodies[2 + 3 * (j + masses.size() * i)] = rr[j * 6 + 2];
    }

    /// for energy fix need change moon center
    std::vector<current_type> moon_res_i(6, 0);
    moonToGeocentric(rr, moon_res_i);

    //------------center mass-------------

    std::vector<current_type> center_mass(3, 0);

    for (int t = 0; t < masses.size(); t++) {
      center_mass[0] += rr[t * 6] * masses[t];
      center_mass[1] += rr[t * 6 + 1] * masses[t];
      center_mass[2] += rr[t * 6 + 2] * masses[t];
    }

    center_mass[0] /= total_mass;
    center_mass[1] /= total_mass;
    center_mass[2] /= total_mass;
    //----------END center mass-------------

    //--------------energy-------------------
    current_type energy = 0;
    current_type k_energy = 0;
    for (int t = 0; t < masses.size(); t++) {
      k_energy +=
          masses[t] *
          (rr[t * 6 + 3] * rr[t * 6 + 3] + rr[t * 6 + 4] * rr[t * 6 + 4] +
           rr[t * 6 + 5] * rr[t * 6 + 5]) /
          2.0;
    }
    current_type p_energy = 0;
    for (int t = 0; t < masses.size(); t++) {
      for (int j = t + 1; j < masses.size(); j++) {
        p_energy +=
            masses[t] * masses[j] /
            std::sqrt((rr[t * 6] - rr[j * 6]) * (rr[t * 6] - rr[j * 6]) +
                       (rr[t * 6 + 1] - rr[j * 6 + 1]) *
                           (rr[t * 6 + 1] - rr[j * 6 + 1]) +
                       (rr[t * 6 + 2] - rr[j * 6 + 2]) *
                           (rr[t * 6 + 2] - rr[j * 6 + 2]));
      }
    }
    energy = k_energy - p_energy;
    //--------END energy---------------------

    //-----------impulse moment---------------
    std::vector<current_type> impulse_moment(3, 0);
    for (int t = 0; t < masses.size(); t++) {
      current_type temp_x = rr[t * 6 + 3] * masses[t];
      current_type temp_y = rr[t * 6 + 4] * masses[t];
      current_type temp_z = rr[t * 6 + 5] * masses[t];

      impulse_moment[0] += rr[t * 6 + 1] * temp_z - rr[t * 6 + 2] * temp_y;
      impulse_moment[1] += rr[t * 6 + 2] * temp_x - rr[t * 6] * temp_z;
      impulse_moment[2] += rr[t * 6] * temp_y - rr[t * 6 + 1] * temp_x;
    }
    //--------------END impulse moment-----------

    data_energy[i] = (abs((energy - init_energy) / init_energy));

    data_impulse_moment[i] =
        std::sqrt((impulse_moment[0] - init_impulse_moment[0]) *
                       (impulse_moment[0] - init_impulse_moment[0]) +
                   (impulse_moment[1] - init_impulse_moment[1]) *
                       (impulse_moment[1] - init_impulse_moment[1]) +
                   (impulse_moment[2] - init_impulse_moment[2]) *
                       (impulse_moment[2] - init_impulse_moment[2])) /
        (init_impulse_moment[0] * init_impulse_moment[0] +
         init_impulse_moment[1] * init_impulse_moment[1] +
         init_impulse_moment[2] * init_impulse_moment[2]);

    data_center[i] = std::sqrt((init_center_mass[0] - center_mass[0]) *
                                    (init_center_mass[0] - center_mass[0]) +
                                (init_center_mass[1] - center_mass[1]) *
                                    (init_center_mass[1] - center_mass[1]) +
                                (init_center_mass[2] - center_mass[2]) *
                                    (init_center_mass[2] - center_mass[2]));

    rr[moonNum * 6] = moon_res_i[0];
    rr[moonNum * 6 + 1] = moon_res_i[1];
    rr[moonNum * 6 + 2] = moon_res_i[2];
    rr[moonNum * 6 + 3] = moon_res_i[3];
    rr[moonNum * 6 + 4] = moon_res_i[4];
    rr[moonNum * 6 + 5] = moon_res_i[5];
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "time: " << duration.count() << " microseconds \n";

  json xyz_data;
  xyz_data["coords"] = data_bodies;
  xyz_data["iter"] = iterations;
  xyz_data["num"] = masses.size();

  writer<current_type> w;
  w.writeRes("log.json", data_energy);
  w.writeRes("bodies.json", xyz_data);
  w.writeRes("moment.json", data_impulse_moment);
  w.writeRes("center.json", data_center);

#ifdef NUMBER_DOUBLE_DOUBLE
  fpu_fix_end(&oldcw);
#endif
  return 0;
}
