#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#define NUMBER_DOUBLE 1
//#define NUMBER_DOUBLE_DOUBLE 1

#include "Integration/methods.h"
#include "Nbodies/energy.h"
#include "Nbodies/summation.h"
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

int main() {
#ifdef NUMBER_DOUBLE_DOUBLE
  unsigned int oldcw;
  fpu_fix_start(&oldcw);
#endif

  using current_type = long double;

  std::ifstream infile("../NBodies/points.txt");
  long double x, y, z, vx, vy, vz, m;
  std::vector<Body<current_type>> bodies;
  std::string name;
  while (infile >> name >> m >> x >> y >> z >> vx >> vy >> vz) {
    bodies.push_back(Body<current_type>({x, y, z}, {vx, vy, vz}, m));
    std::cout << "parsed " << name << " : " << x << " " << y << " " << z << " "
              << vx << " " << vy << " " << vz << " " << m << std::endl;
  }
  std::cout << bodies[0].r.X << std::endl;
  std::cout << "hello";

  /// for energy fix need change moon center
  current_type earth_x = bodies[earthNum].r.X, earth_y = bodies[earthNum].r.Y,
               earth_z = bodies[earthNum].r.Z;
  current_type earth_vx = bodies[earthNum].v.X, earth_vy = bodies[earthNum].v.Y,
               earth_vz = bodies[earthNum].v.Z;
  current_type gcmoon_x = bodies[moonNum].r.X, gcmoon_y = bodies[moonNum].r.Y,
               gcmoon_z = bodies[moonNum].r.Z;
  current_type gcmoon_vx = bodies[moonNum].v.X, gcmoon_vy = bodies[moonNum].v.Y,
               gcmoon_vz = bodies[moonNum].v.Z;

  bodies[moonNum].r.X = earth_x + gcmoon_x;
  bodies[moonNum].r.Y = earth_y + gcmoon_y;
  bodies[moonNum].r.Z = earth_z + gcmoon_z;
  bodies[moonNum].v.X = earth_vx + gcmoon_vx;
  bodies[moonNum].v.Y = earth_vy + gcmoon_vy;
  bodies[moonNum].v.Z = earth_vz + gcmoon_vz;

  current_type total_mass =
      summation<current_type, total_mass_proxy<current_type>>(
          total_mass_proxy<current_type>(bodies), bodies.size());
  vec<current_type> init_center_mass =
      summation<vec<current_type>,
                mass_center_proxy<vec<current_type>, current_type>>(
          mass_center_proxy<vec<current_type>, current_type>(bodies),
          bodies.size());
  init_center_mass = init_center_mass / total_mass;
  vec<current_type> init_vel_mass =
      summation<vec<current_type>,
                mass_vel_proxy<vec<current_type>, current_type>>(
          mass_vel_proxy<vec<current_type>, current_type>(bodies),
          bodies.size());
  init_vel_mass = init_vel_mass / total_mass;

  std::cout << init_center_mass.X << " " << init_center_mass.Y << " "
            << init_center_mass.Z << " " << init_vel_mass << std::endl;

  for (int i = 0; i < bodies.size(); i++) {
    // bodies[i].r -= init_center_mass;
    // bodies[i].v -= init_vel_mass;
  }

  current_type init_energy =
      summation<current_type, kinetic_energy_proxy<current_type>>(
          kinetic_energy_proxy(bodies), bodies.size()) /
          current_type(2) +
      summation<current_type, potential_energy_proxy<current_type>>(
          potential_energy_proxy(bodies), bodies.size() * bodies.size()) /
          current_type(2);

  vec<current_type> init_impulse_moment =
      summation<vec<current_type>,
                impulse_moment_proxy<vec<current_type>, current_type>>(
          impulse_moment_proxy<vec<current_type>, current_type>(bodies),
          bodies.size());

  bodies[moonNum].r.X = gcmoon_x;
  bodies[moonNum].r.Y = gcmoon_y;
  bodies[moonNum].r.Z = gcmoon_z;
  bodies[moonNum].v.X = gcmoon_vx;
  bodies[moonNum].v.Y = gcmoon_vy;
  bodies[moonNum].v.Z = gcmoon_vz;

  std::vector<current_type> data_bodies, data_energy, data_impulse_moment,
      data_center;

  current_type h(0.0625);
  int iterations = 5840;

  data_bodies.resize(bodies.size() * iterations * 3);
  data_energy.resize(iterations);
  data_impulse_moment.resize(iterations);
  data_center.resize(iterations);

  std::vector<current_type> coefs = initDDCoef<current_type>();
  auto start = std::chrono::high_resolution_clock::now();

  std::vector<current_type> masses;
  masses.resize(bodies.size());
  for (int i = 0; i < bodies.size(); i++) {
    masses[i] = bodies[i].m;
  }
  for (int i = 0; i < iterations; i++) {
    std::vector<current_type> x;
    x.resize(6 * bodies.size());
    for (int j = 0; j < bodies.size(); j++) {
      x[j * 6] = bodies[j].r.X;
      x[j * 6 + 1] = bodies[j].r.Y;
      x[j * 6 + 2] = bodies[j].r.Z;
      x[j * 6 + 3] = bodies[j].v.X;
      x[j * 6 + 4] = bodies[j].v.Y;
      x[j * 6 + 5] = bodies[j].v.Z;
    }
    Euler(x, masses, h);

    for (int j = 0; j < bodies.size(); j++) {
      bodies[j].r.X = x[j * 6];
      bodies[j].r.Y = x[j * 6 + 1];
      bodies[j].r.Z = x[j * 6 + 2];
      bodies[j].v.X = x[j * 6 + 3];
      bodies[j].v.Y = x[j * 6 + 4];
      bodies[j].v.Z = x[j * 6 + 5];
    }

    for (int j = 0; j < bodies.size(); j++) {
      data_bodies[3 * (j + bodies.size() * i)] = bodies[j].r.X;
      data_bodies[1 + 3 * (j + bodies.size() * i)] = bodies[j].r.Y;
      data_bodies[2 + 3 * (j + bodies.size() * i)] = bodies[j].r.Z;
    }

    /// for energy fix need change moon center
    current_type earth_x = bodies[earthNum].r.X, earth_y = bodies[earthNum].r.Y,
                 earth_z = bodies[earthNum].r.Z;
    current_type earth_vx = bodies[earthNum].v.X,
                 earth_vy = bodies[earthNum].v.Y,
                 earth_vz = bodies[earthNum].v.Z;
    current_type gcmoon_x = bodies[moonNum].r.X, gcmoon_y = bodies[moonNum].r.Y,
                 gcmoon_z = bodies[moonNum].r.Z;
    current_type gcmoon_vx = bodies[moonNum].v.X,
                 gcmoon_vy = bodies[moonNum].v.Y,
                 gcmoon_vz = bodies[moonNum].v.Z;

    bodies[moonNum].r.X = earth_x + gcmoon_x;
    bodies[moonNum].r.Y = earth_y + gcmoon_y;
    bodies[moonNum].r.Z = earth_z + gcmoon_z;
    bodies[moonNum].v.X = earth_vx + gcmoon_vx;
    bodies[moonNum].v.Y = earth_vy + gcmoon_vy;
    bodies[moonNum].v.Z = earth_vz + gcmoon_vz;

    vec<current_type> center_mass =
        summation<vec<current_type>,
                  mass_center_proxy<vec<current_type>, current_type>>(
            mass_center_proxy<vec<current_type>, current_type>(bodies),
            bodies.size());
    center_mass = center_mass / total_mass;

    current_type energy =
        summation<current_type, kinetic_energy_proxy<current_type>>(
            kinetic_energy_proxy(bodies), bodies.size()) /
            current_type(2) +
        summation<current_type, potential_energy_proxy<current_type>>(
            potential_energy_proxy(bodies), bodies.size() * bodies.size()) /
            current_type(2);

    vec<current_type> impulse_moment =
        summation<vec<current_type>,
                  impulse_moment_proxy<vec<current_type>, current_type>>(
            impulse_moment_proxy<vec<current_type>, current_type>(bodies),
            bodies.size());

    data_energy[i] = (abs((energy - init_energy) / init_energy));

    data_impulse_moment[i] = (abs((impulse_moment - init_impulse_moment).Len() /
                                  init_impulse_moment.Len()));

    data_center[i] = (abs((init_center_mass - center_mass).Len()));

    bodies[moonNum].r.X = gcmoon_x;
    bodies[moonNum].r.Y = gcmoon_y;
    bodies[moonNum].r.Z = gcmoon_z;
    bodies[moonNum].v.X = gcmoon_vx;
    bodies[moonNum].v.Y = gcmoon_vy;
    bodies[moonNum].v.Z = gcmoon_vz;
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "time: " << duration.count() << " microseconds \n";

  json xyz_data;
  xyz_data["coords"] = data_bodies;
  xyz_data["iter"] = iterations;
  xyz_data["num"] = bodies.size();

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
