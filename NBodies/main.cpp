#include <iostream>

#include "Integration/integration.h"
#include <vector>
#include "Writer/writer.h"
using namespace std;
int main() {

    vector< Body<double> > bodies;


    bodies.push_back({{ 0, 1, 0 }, { 35, 0, 0 }, 1, {0,0,3}});
    bodies.push_back({{ 0, 2, 1 }, { 3,0,0 }, 1, {0,3,1}});
    bodies.push_back({{ 1, 2, 0 }, { 2,0,0 }, 1, {5,2,4}});
    bodies.push_back({{ 4,0,0 }, { 0,13,0 }, 1, {3,6,8}});
    bodies.push_back({{ 0, 2, 0 }, { 9,0,0 }, 1, {1,0,0}});
    bodies.push_back({{ 0, 8, 2 }, { 6, 0, 0 }, 1, {2,1,1}});
    bodies.push_back({{ 1, 1, 1 }, { 54,0,0 }, 1, {1,1,1}});
    bodies.push_back({{ 0, 14, 18 }, { 4,0,0 }, 1, {1,0,1}});

    Integration<double> integrator;
    Utils<double> util;
    vector<double> energy;
    vector<double> iterations;

    for(int i = 0; i < 10000; i++){
        integrator.compute(bodies, 0.001);
        energy.push_back(util.energy(bodies));
        printf("Energy: %.18le \n", energy.back());
        iterations.push_back(i);
    }

    writer<double> w;
    w.writeRes("../log.txt", iterations, energy);

    return 0;
}
