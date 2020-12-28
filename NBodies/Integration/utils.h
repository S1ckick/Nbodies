//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_UTILS_H
#define NBODIES_UTILS_H
#include "../Bodies/bodies.h"

template <typename Type>
class Utils{
public:
    vec<Type> partial_step(vec<Type> f, vec<Type> df, Type scale, double t);
    Type energy(std::vector< Body <Type>> &bodies);
};

template <typename Type>
vec<Type> Utils<Type>::partial_step(vec<Type> f, vec<Type> df, Type scale, double t){
    return f + df * t * scale;
}

template <typename Type>
Type Utils<Type>::energy(std::vector< Body <Type>> &bodies){
    Type KE = 0.0;
    for(int i = 0; i < bodies.size(); i++){
        KE+= bodies[i].m * bodies[i].v.Len2() / 2.0;
    }
    Type PE = 0;
    for(int i = 0; i < bodies.size(); i++){
        for(int j = i + 1; j < bodies.size(); j++){
            vec<Type> dist = bodies[i].r - bodies[j].r;
            PE+=G*bodies[i].m * bodies[j].m / dist.Len();
        }
    }
    return PE+KE;
}

#endif //NBODIES_UTILS_H
