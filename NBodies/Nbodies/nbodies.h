//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_NBODIES_H
#define NBODIES_NBODIES_H

#include "../Utils/vec.h"

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


#endif //NBODIES_NBODIES_H
