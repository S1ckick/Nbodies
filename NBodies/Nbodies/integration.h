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

#include "Nbodies.h"

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
