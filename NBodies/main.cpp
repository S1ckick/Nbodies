#include <iostream>
#include <vector>

#include "Nbodies/integration.h"
#include "Integration/methods.h"
#include "Writer/writer.h"

int main() {

    std::vector<Body<double>> bodies;

    //bodies.push_back(Body<double>({0,3e4,0},{0,0,0},2e14));
    //bodies.push_back(Body<double>({0,3e4 + 1.5e3,0},{3,0,0},6));

    bodies.push_back(Body<double>({0,0,0},{0,0,0},2e14));
    bodies.push_back(Body<double>({ 0, 1.4e3, 0 },{ 3,0,0 },6));
    bodies.push_back(Body<double>({ 0, 1.3e3, 0 },{ 3,0,0 },6));
    bodies.push_back(Body<double>({ 0, 1.5e3, 0 },{ 2,0,sqrt(5) },6));
    bodies.push_back(Body<double>({ 0, 1.2e3, 0 },{ 0,0,3 },6));
    bodies.push_back(Body<double>({ 0, 1.5e3, 0 },{ -2,0,sqrt(5) },6));


    //bodies.push_back(Body<double>({3e3,0,0},{0,0,0},2e14));
    //bodies.push_back(Body<double>({ 4.5e3, 0, 0 },{ 3,0,0 },6));

    json data_energy, data_bodies;

    for (int i = 0; i < 100000; i++) {
        RungeKutta4(bodies, 0.1);

        if( true ) {
            for(int j = 0; j < bodies.size(); j++){
                data_bodies[j]["X"][i] = bodies[j].r.X;
                data_bodies[j]["Y"][i] = bodies[j].r.Y;
                data_bodies[j]["Z"][i] = bodies[j].r.Z;
            }

            data_energy["n"].push_back(i);
            data_energy["energy"].push_back(Energy(bodies));

            printf("Energy: %.64le \n", (double)data_energy["energy"][i]);
        }
    }

    writer<double> w;
    w.writeRes("../log.json", data_energy);
    w.writeRes("../bodies.json", data_bodies);


    return 0;
}
