#include <iostream>
#include <vector>

#include "Nbodies/summation.h"
#include "Integration/methods.h"
#include "Writer/writer.h"


#include <qd/dd_real.h>
#include <qd/fpu.h>
#include <string>

using namespace std;


#include <chrono>
int main() {

    unsigned int oldcw;
    fpu_fix_start(&oldcw);

    using current_type = dd_real;

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Body<current_type>> bodies;

    //bodies.push_back(Body<double>({0,3e4,0},{0,0,0},2e14));
    //bodies.push_back(Body<double>({0,3e4 + 1.5e3,0},{3,0,0},6));

   //make_universe(bodies,5,current_type(100),current_type(100),current_type(100));


    bodies.push_back(Body<current_type>({current_type(0),current_type(0),current_type(0)},{current_type(0),current_type(0),current_type(0)},current_type(2e14)));
    bodies.push_back(Body<current_type>({ current_type(0), current_type(1.4e3), current_type(0) },{ current_type(3),current_type(0),current_type(0) },current_type(6)));
    bodies.push_back(Body<current_type>({ current_type(0), current_type(1.3e3), current_type(0) },{ current_type(3),current_type(0),current_type(0) },current_type(6)));
    bodies.push_back(Body<current_type>({ current_type(0), current_type(1.2e3), current_type(0) },{ current_type(3),current_type(0),current_type(0) },current_type(6)));
    bodies.push_back(Body<current_type>({ current_type(0), current_type(1.1e3), current_type(0) },{ current_type(3),current_type(0),current_type(0) },current_type(6)));
    bodies.push_back(Body<current_type>({ current_type(0), current_type(1.56e3), current_type(0) },{ current_type(3),current_type(0),current_type(0) },current_type(6)));
    bodies.push_back(Body<current_type>({ current_type(0), current_type(1.5e3), current_type(0) },{ current_type(3),current_type(0),current_type(0) },current_type(6)));



    //bodies.push_back(Body<double>({3e3,0,0},{0,0,0},2e14));
    //bodies.push_back(Body<double>({ 4.5e3, 0, 0 },{ 3,0,0 },6));

    json data_energy, data_bodies;

    current_type init_energy = summation<current_type, kinetic_energy_proxy<current_type>>(kinetic_energy_proxy(bodies), bodies.size()) / current_type(2)+
                         summation<current_type, potential_energy_proxy<current_type>>(potential_energy_proxy(bodies), bodies.size() * bodies.size()) / current_type(2);

    for (int i = 0; i < 100000; i++) {
         dormanPrince8(bodies, current_type(0.1));

        if( true ) {
            for(int j = 0; j < bodies.size(); j++){

                    data_bodies[j]["X"][i] = bodies[j].r.X.to_string(32, 0, 0, false, false);
                    data_bodies[j]["Y"][i] = bodies[j].r.Y.to_string(32, 0, 0, false, false);
                    data_bodies[j]["Z"][i] = bodies[j].r.Z.to_string(32, 0, 0, false, false);


                    /*
                    data_bodies[j]["X"][i] = bodies[j].r.X;
                    data_bodies[j]["Y"][i] = bodies[j].r.Y;
                    data_bodies[j]["Z"][i] = bodies[j].r.Z;
                    */

            }

            data_energy["n"].push_back(i);
            current_type energy = summation<current_type, kinetic_energy_proxy<current_type>>(kinetic_energy_proxy(bodies), bodies.size()) / current_type(2) +
                            summation<current_type, potential_energy_proxy<current_type>>(potential_energy_proxy(bodies),bodies.size() * bodies.size()) / current_type(2);

            data_energy["energy"].push_back(abs((energy - init_energy)/init_energy).to_string(32, 0, 0, false, false));
            //data_energy["energy"].push_back(abs((energy - init_energy)/init_energy));

            //printf("Energy: %.64le \n", (double)data_energy["energy"][i]);
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(stop - start);
    printf("time: %lld microseconds \n", duration.count());
    writer<current_type> w;
    w.writeRes("../log.json", data_energy);
    w.writeRes("../bodies.json", data_bodies);

    fpu_fix_end(&oldcw);
    return 0;
}
