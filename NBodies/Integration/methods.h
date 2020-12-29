//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_METHODS_H
#define NBODIES_METHODS_H

#include "../Nbodies/nbodies.h"

template<typename Type>
void copyBodies(const std::vector<Body<Type>> &bodies, std::vector<Body<Type>> &newBodies){
    newBodies.clear();
    for(int i = 0; i < bodies.size(); i++){
        newBodies.push_back(bodies[i]);
    }
}

template <typename Type>
void f(const std::vector<Body<Type>> &bodies, std::vector<Body<Type>> &fbodies){
    for(int body_1_idx = 0; body_1_idx < bodies.size(); body_1_idx++){
        vec<Type> total_force(Type(0),Type(0),Type(0));
        for(int body_2_idx = 0; body_2_idx < bodies.size(); body_2_idx++){
            if(body_1_idx == body_2_idx) continue;

            total_force += bodies[body_1_idx].IteractSubtotalForce(bodies[body_2_idx]);
        }
        fbodies[body_1_idx].r = bodies[body_1_idx].v;
        fbodies[body_1_idx].v = total_force;
    }
}


template <typename Type>
std::vector<Body<Type>> addBodies(const std::vector<Body<Type>> &bodies, const std::vector<Body<Type>> &fbodies){
    std::vector<Body<Type>> newBodies;
    copyBodies(bodies, newBodies);
    for(int i = 0; i < bodies.size(); i++){
        newBodies[i].r = bodies[i].r + fbodies[i].r;
        newBodies[i].v = bodies[i].v + fbodies[i].v;
    }
    return newBodies;
}

template <typename Type>
std::vector<Body<Type>> multBodies(const std::vector<Body<Type>> &bodies, Type number){
    std::vector<Body<Type>> newBodies;
    copyBodies(bodies, newBodies);

    for(int i = 0; i < bodies.size(); i++){
        newBodies[i].r = bodies[i].r * number;
        newBodies[i].v = bodies[i].v * number;
    }
    return newBodies;
}



template <typename Type>
vec<Type> RungeKutta4(std::vector<Body<Type>> &bodies, double h) {

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
