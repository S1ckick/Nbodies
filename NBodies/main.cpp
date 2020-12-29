#include <iostream>
#include <vector>

#include "Nbodies/integration.h"
#include "Integration/methods.h"
#include "Writer/writer.h"

int main() {

    std::vector<Body<double>> b;


    b.push_back(Body<double>({0,0,0},{0,0,0},2e14));
    b.push_back(Body<double>({ 0, 1.5e3, 0 },{ 3,0,0 },6));

    std::vector<double> energy;
    std::vector<double> iterations;
    std::vector<double> x_1;
    std::vector<double> x_2;
    std::vector<double> y_1;
    std::vector<double> y_2;
    std::vector<double> z_1;
    std::vector<double> z_2;



    for (int i = 0; i < 100000; i++) {
        RungeKutta4(b, 0.1);

        if(true) {
            energy.push_back(Energy(b));
            iterations.push_back(i);
            x_1.push_back(b[0].r.X);
            y_1.push_back(b[0].r.Y);
            z_1.push_back(b[0].r.Z);

            x_2.push_back(b[1].r.X);
            y_2.push_back(b[1].r.Y);
            z_2.push_back(b[1].r.Z);

            printf("Energy: %.64le \n", energy.back());
        }
    }

    writer<double> w;
    w.writeRes("../log.txt", iterations, energy);
    w.writeBodies("../bodies.txt", x_1, y_1, z_1, x_2, y_2, z_2);


    return 0;
}
