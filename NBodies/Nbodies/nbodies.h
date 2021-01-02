//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_NBODIES_H
#define NBODIES_NBODIES_H

#include "../Utils/vec.h"
#include "summation.h"

static const double Gamma = 6.6743e-11;

template <class Type>
struct Body {
    Type m;
    vec<Type> r, v;

    Body(vec<Type> r1, vec<Type> v1, Type m1)
            : r(r1), v(v1), m(m1) {};

    vec<Type> IteractSubtotalForce(const Body<Type> &b2) const{
        vec<Type> dist = r - b2.r;
        Type len = dist.Len();
        if (len < 1e-8){
            len = 1e-8;
        }
        return dist * (-(Gamma * b2.m / (len * len * len)));
    }

    Body& operator=(Body body){
        m = body.m;
        r = body.r;
        v = body.v;
    }

    Body& operator+(Body body){
        r = body.r;
        v = body.v;
    }

    Body& operator*(Type number){
        r = r * number;
        v = v * number;
    }


};

template<typename Type>
void copyBodies(const std::vector<Body<Type>> &bodies, std::vector<Body<Type>> &newBodies){
    newBodies.clear();
    for(int i = 0; i < bodies.size(); i++){
        newBodies.push_back(bodies[i]);
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
std::vector<Body<Type>> multBodies(const std::vector<Body<Type>> &bodies, double number){
    std::vector<Body<Type>> newBodies;
    copyBodies(bodies, newBodies);

    for(int i = 0; i < bodies.size(); i++){
        newBodies[i].r = bodies[i].r * number;
        newBodies[i].v = bodies[i].v * number;
    }
    return newBodies;
}

template <typename Type>
struct acceleration_proxy {
  const std::vector<Body<Type>> &bodies;
  const size_t bodyind;
  explicit acceleration_proxy(const std::vector<Body<Type>> &bodies, size_t bi)
      : bodies(bodies), bodyind(bi) {}
  vec<Type> operator[](size_t n) const { 
    if (n == bodyind)
      return vec<Type>(0);
    return bodies[bodyind].IteractSubtotalForce(bodies[n]); 
  }
};

// derivative
template <typename Type>
void f(const std::vector<Body<Type>> &bodies,
       std::vector<Body<Type>> &fbodies) {
  for (int body_1_idx = 0; body_1_idx < bodies.size(); body_1_idx++) {
    acceleration_proxy<Type> a(bodies, body_1_idx);
    vec<Type> total_force = 
      summation<vec<Type>, acceleration_proxy<Type>>(a, bodies.size());
    //vec<Type> total_force(Type(0), Type(0), Type(0));
    //for (int body_2_idx = 0; body_2_idx < bodies.size(); body_2_idx++) {
    //  if (body_1_idx == body_2_idx) continue;
    //
    //  total_force +=
    //      bodies[body_1_idx].IteractSubtotalForce(bodies[body_2_idx]);
    //}
    fbodies[body_1_idx].r = bodies[body_1_idx].v;
    fbodies[body_1_idx].v = total_force;
  }
}

#endif //NBODIES_NBODIES_H
