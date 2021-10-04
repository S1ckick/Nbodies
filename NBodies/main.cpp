#include <fstream>
#include <iostream>
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
#endif

#include <string>

#include "Utils/helper.h"

using namespace std;

#include <chrono>
int main() {
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
    bodies.push_back(Body<current_type>({x, y, z}, {vx, vy, vz}, m));
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

  for (int i = 0; i < bodies.size(); i++) {
    bodies[i].r -= init_center_mass;
    bodies[i].v -= init_vel_mass;
  }

  json data_energy, data_impulse_moment, data_center, data_bodies;
  current_type h(0.3);
  int iterations = 100000;

  std::vector<current_type> coefs = initDDCoef<current_type>();
  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < iterations; i++) {
    RungeKutta4(bodies, h);

    for (int j = 0; j < bodies.size(); j++) {
#ifdef NUMBER_DOUBLE_DOUBLE
      data_bodies[j]["X"][i] = bodies[j].r.X._hi();
      data_bodies[j]["Y"][i] = bodies[j].r.Y._hi();
      data_bodies[j]["Z"][i] = bodies[j].r.Z._hi();
#endif
#ifdef NUMBER_DOUBLE

      data_bodies[j]["X"][i] = bodies[j].r.X;
      data_bodies[j]["Y"][i] = bodies[j].r.Y;
      data_bodies[j]["Z"][i] = bodies[j].r.Z;
#endif
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

    data_energy["n"].push_back(i);
#ifdef NUMBER_DOUBLE_DOUBLE
    data_energy["energy"].push_back(
        (abs((energy - init_energy) / init_energy)).to_string());
#endif
#ifdef NUMBER_DOUBLE
    data_energy["energy"].push_back(abs((energy - init_energy) / init_energy));
#endif

    data_impulse_moment["n"].push_back(i);
#ifdef NUMBER_DOUBLE_DOUBLE
    data_impulse_moment["moment"].push_back(
        abs((impulse_moment - init_impulse_moment).Len() /
            init_impulse_moment.Len())
            .to_string());
#endif
#ifdef NUMBER_DOUBLE
    data_impulse_moment["moment"].push_back(
        abs((impulse_moment - init_impulse_moment).Len() /
            init_impulse_moment.Len()));
#endif

    data_center["n"].push_back(i);
#ifdef NUMBER_DOUBLE_DOUBLE
    data_center["center"].push_back(
        abs((init_center_mass - center_mass).Len()).to_string());
#endif
#ifdef NUMBER_DOUBLE
    data_center["center"].push_back(
        abs((init_center_mass - center_mass).Len()));
#endif
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = duration_cast<std::chrono::microseconds>(stop - start);
  printf("time: %lld microseconds \n", duration.count());

  writer<current_type> w;
  w.writeRes("../log.json", data_energy);
  w.writeRes("../bodies.json", data_bodies);
  w.writeRes("../moment.json", data_impulse_moment);
  w.writeRes("../center.json", data_center);

#ifdef NUMBER_DOUBLE_DOUBLE
  fpu_fix_end(&oldcw);
#endif
  return 0;
}
