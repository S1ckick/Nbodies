#include <iostream>

#include "Integration/integration.h"
#include <vector>
using namespace std;
int main() {

    vector< Body<double> > bodies;

    bodies.push_back({{ 0, 0, 0 }, { 0,0,0 }, 2e30});
    bodies.push_back({{ 0, 5.0e10, 0 }, { 4, 0, 0 }, 3.285e23});
    bodies.push_back({{ 0, 1.1e11, 0 }, { 35, 0, 0 }, 4.8e24});
    bodies.push_back({{ 0, 1.5e11, 0 }, { 3,0,0 }, 6e24});
    bodies.push_back({{ 0, 2.2e11, 0 }, { 2,0,0 }, 2.4e24});
    bodies.push_back({{ 7.7e11,0,0 }, { 0,13,0 }, 1e28});
    bodies.push_back({{ 0, 1.4e12, 0 }, { 9,0,0 }, 5.7e26});
    bodies.push_back({{ 0, 2.8e12, 0 }, { 6, 0, 0 }, 8.7e25});
    bodies.push_back({{ 0, 4.5e12, 0 }, { 54,0,0 }, 1e26});
    bodies.push_back({{ 0, 7.3e12, 0 }, { 4,0,0 }, 1.3e22});

    Integration<double> integrator;
    for(int i = 0; i < 10; i++){
        integrator.compute(bodies, 1);
        for(int j = 0; j< bodies.size(); j++){
            cout << j << ") pos: " <<bodies[j].r.X << " " << bodies[j].r.Y << " " << bodies[j].r.Z << endl;
            cout <<"  vel: "<<bodies[j].v.X << " " << bodies[j].v.Y << " " << bodies[j].v.Z << endl;
            cout <<"  acc: "<<bodies[j].a.X << " " << bodies[j].a.Y << " " << bodies[j].a.Z << endl;
            cout << "------------------------------- \n";
        }
    }

    return 0;
}
