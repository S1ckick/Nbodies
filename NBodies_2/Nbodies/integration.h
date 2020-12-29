//
// Created by Максим on 29.12.2020.
//

#ifndef NBODIES_INTEGRATION_H
#define NBODIES_INTEGRATION_H

#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <string>

#include "nbodies.h"

template <typename Type>
void Compute(std::vector<Body<Type>> &bodies, Type step,
             vec<Type> (*Method)(const Body<Type> &, double)) {
    for (size_t i = 0; i < bodies.size(); i++) {
        Type mass1 = bodies[i].m;
        vec<Type> acceleration(Type(0),Type(0),Type(0));
        for (size_t j = 0; j < bodies.size(); j++) {
            if (i == j) continue;

            Type mass2 = bodies[j].m, da = bodies[i].IteractSubtotalForce(bodies[j]);
            vec<Type> dist = bodies[j].r - bodies[i].r;
            Body<Type> tmp(dist, bodies[i].v, 0, vec<Type>(da));
            acceleration += Method(tmp, step);
        }
        bodies[i].Update(step, acceleration);
    }
}

template <typename Type>
Type Energy(std::vector<Body<Type>> &bodies) {
    Type KE = 0.0;
    for (int i = 0; i < bodies.size(); i++) {
        KE += bodies[i].m * bodies[i].v.Len2() / 2.0;
    }
    Type PE = 0;
    for (int i = 0; i < bodies.size(); i++) {
        for (int j = 0; j < bodies.size(); j++) {
            if (i == j) continue;
            vec<Type> dist = bodies[i].r - bodies[j].r;
            PE += Gamma * bodies[i].m * bodies[j].m / dist.Len();
        }
    }
    return PE + KE;
}



#endif //NBODIES_INTEGRATION_H
