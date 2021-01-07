
#ifndef NBODIES_ENERGY_H
#define NBODIES_ENERGY_H

#ifdef NUMBER_DOUBLE_DOUBLE
#include <qd/dd_real.h>
#endif

#include <vector>
#include <string>

#include "nbodies.h"

template<typename Type>
Type potential_energy(const std::vector<Body<Type>> &bodies, size_t body1, size_t body2){
    vec<Type> dr = bodies[body1].r - bodies[body2].r;
    Type r2 = dr.Len2();
    return -(Gamma<Type> * bodies[body1].m * bodies[body2].m ) / sqrt(r2);
}

template<typename Type>
struct potential_energy_proxy{
    const std::vector<Body<Type>> &bodies;
    explicit potential_energy_proxy(const std::vector<Body<Type>> &bodies) : bodies(bodies) {};
    Type operator [](size_t n) const{
        size_t n1 = n / bodies.size();
        size_t n2 = n % bodies.size();
        if (n1 == n2){
            return Type(0);
        }
        return potential_energy<Type>(bodies, n1, n2);
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


template <class Type, typename BodyType>
struct impulse_moment_proxy{
    const std::vector<Body<BodyType>> &bodies;
    explicit impulse_moment_proxy(const std::vector<Body<BodyType>> &bodies) : bodies(bodies) {};
    Type operator[](size_t n) const {
        return bodies[n].r % (bodies[n].v * bodies[n].m);
    }
};

template <class Type, typename BodyType>
struct mass_center_proxy {
    const std::vector<Body<BodyType>> &bodies;
    explicit mass_center_proxy(const std::vector<Body<BodyType>> &bodies) : bodies(bodies) {};
    Type operator[](size_t n) const {
        return bodies[n].r * bodies[n].m;
    }
};

template <class Type, typename BodyType>
struct mass_vel_proxy {
    const std::vector<Body<BodyType>> &bodies;
    explicit mass_vel_proxy(const std::vector<Body<BodyType>> &bodies) : bodies(bodies) {};
    Type operator[](size_t n) const {
        return bodies[n].v * bodies[n].m;
    }
};

#endif //NBODIES_ENERGY_H
