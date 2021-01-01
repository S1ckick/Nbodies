//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_METHODS_H
#define NBODIES_METHODS_H

#include "../Nbodies/nbodies.h"

template <typename Type>
void RungeKutta4(std::vector<Body<Type>> &bodies, double h) {

    std::vector<Body<Type>> k_1;
    copyBodies(bodies, k_1);

    std::vector<Body<Type>> temp;
    copyBodies(bodies, temp);


    f(temp,k_1);

    //temp = bodies + 0.5*k1*h
    copyBodies(addBodies(multBodies(multBodies(k_1, 0.5), h), temp), temp);

    std::vector<Body<Type>> k_2;
    copyBodies(bodies, k_2);
    f(temp, k_2);

    //temp += 0.5*k2*h
    copyBodies(addBodies(multBodies(multBodies(k_2, 0.5), h),temp), temp);

    std::vector<Body<Type>> k_3;
    copyBodies(bodies, k_3);
    f(temp, k_3);

    //temp +=  1.0*k3*h
    copyBodies(addBodies(multBodies(multBodies(k_3, 1.0), h),temp), temp);

    std::vector<Body<Type>> k_4;
    copyBodies(bodies, k_4);
    f(temp, k_4);

    copyBodies(multBodies(k_1, 1.0/6.0), k_1);
    copyBodies(multBodies(k_2, 1.0/3.0), k_2);
    copyBodies(multBodies(k_3, 1.0/3.0), k_3);
    copyBodies(multBodies(k_4, 1.0/8.0), k_4);

    //y += 	dt( k_1/6 + k_2/3 + k_3/3 + k_4/8 )
    copyBodies(addBodies(multBodies(addBodies(addBodies(addBodies(k_1, k_2), k_3), k_4), h), bodies), bodies);
}




#endif //NBODIES_METHODS_H
