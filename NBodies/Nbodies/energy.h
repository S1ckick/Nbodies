//
// Created by Максим on 29.12.2020.
//

#ifndef NBODIES_ENERGY_H
#define NBODIES_ENERGY_H

#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <string>

#include "nbodies.h"

template<typename Type>
Type potential_energy(const std::vector<Body<Type>> &bodies, size_t body1, size_t body2){
    vec<Type> dr = bodies[body1].r - bodies[body2].r;
    Type r2 = dr.Len2();
    if(r2 < 1e-8){
        return Type(0);
    }
    return - ( Gamma * bodies[body1].m * bodies[body2].m ) / sqrt(r2);
}

template<typename Type>
struct potential_energy_proxy{
    const std::vector<Body<Type>> &bodies;
        explicit potential_energy_proxy(const std::vector<Body<Type>> &bodies) : bodies(bodies) {};
        Type operator [](size_t n) const{
            size_t n1 = n / bodies.size();
            size_t n2 = n % bodies.size();
            if(n1 == n2){
                return Type(0);
            }
            return potential_energy(bodies, n1, n2);
        }

    };

template <typename Type>
struct kinetic_energy_proxy{
    const std::vector<Body<Type>> &bodies;
    explicit kinetic_energy_proxy(const std::vector<Body<Type>> &bodies) : bodies(bodies) {};
    Type operator[](size_t n) const{
        return bodies[n].v.Len2() * bodies[n].m;
    }
};



#endif //NBODIES_ENERGY_H
