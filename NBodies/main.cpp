#include <iostream>

#include "Integration/integration.h"
#include <vector>
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
    for(int i = 0; i < 10; i++){

        for(int j = 0; j< bodies.size(); j++){
            cout << j << ") pos: " <<bodies[j].r.X << " " << bodies[j].r.Y << " " << bodies[j].r.Z << endl;
            cout <<"  vel: "<<bodies[j].v.X << " " << bodies[j].v.Y << " " << bodies[j].v.Z << endl;
            cout <<"  acc: "<<bodies[j].a.X << " " << bodies[j].a.Y << " " << bodies[j].a.Z << endl;
            cout << "------------------------------- \n";
        }
        integrator.compute(bodies, 0.001);
        double energy = util.energy(bodies);
        printf("Energy: %.18le \n", energy);
    }

    return 0;
}
