#include <fstream>
#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <filesystem>

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
  std::cout << std::filesystem::current_path();
#ifdef NUMBER_DOUBLE_DOUBLE
  unsigned int oldcw;
  fpu_fix_start(&oldcw);
#endif

  using current_type = long double;

  std::ifstream infile("../points.txt");
  long double x, y, z, vx, vy, vz, m;
  std::vector<Body<current_type>> bodies;
  std::string name;
  while (infile >> name >> m >> x >> y >> z >> vx >> vy >> vz) {
    bodies.push_back(Body<current_type>({vx, vy, vz}, {x, y, z}, m));
    std::cout << "parsed " << name << " : " << x << " " << y << " " << z << " "
              << vx << " " << vy << " " << vz << " " << m << std::endl;
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

  std::vector<current_type> data_bodies, data_energy, data_impulse_moment, data_center;

  current_type h(0.3);
  int iterations = 100000;

  data_bodies.resize(bodies.size() * iterations * 3);
  data_energy.resize(iterations);
  data_impulse_moment.resize(iterations);
  data_center.resize(iterations);

  std::vector<current_type> coefs = initDDCoef<current_type>();
  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < iterations; i++) {
    RungeKutta4(bodies, h);

    for (int j = 0; j < bodies.size(); j++) {
      data_bodies[3 * (j + bodies.size() * i)] = bodies[j].r.X;
      data_bodies[1 + 3 * (j + bodies.size() * i)] = bodies[j].r.Y;
      data_bodies[2 + 3 * (j + bodies.size() * i)] = bodies[j].r.Z;
    }

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

    data_impulse_moment[i] = (
        abs((impulse_moment - init_impulse_moment).Len() /
            init_impulse_moment.Len()));

    data_center[i] = (abs((init_center_mass - center_mass).Len()));
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "time: " << duration.count() << " microseconds \n";

  json xyz_data;
  xyz_data["coords"] = data_bodies;
  xyz_data["iter"] = iterations;
  xyz_data["num"] = bodies.size();

  writer<current_type> w;
  w.writeRes("../log.json", data_energy);
  w.writeRes("../bodies.json", xyz_data);
  w.writeRes("../moment.json", data_impulse_moment);
  w.writeRes("../center.json", data_center);

#ifdef NUMBER_DOUBLE_DOUBLE
  fpu_fix_end(&oldcw);
#endif
  return 0;
}
