//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_NBODIES_H
#define NBODIES_NBODIES_H

#include "../Utils/vec.h"

const double Gamma = 6.6743e-11;

template <class Type>
struct Body {
    Type m;
    vec<Type> r, v, a;

    Body(vec<Type> r1, vec<Type> v1, Type m1, vec<Type> a1 = {Type(0), Type(0), Type(0)})
            : r(r1), v(v1), a(a1), m(m1) {};


    vec<Type> IteractSubtotalForce(const Body<Type> &b2) const{
        vec<Type> dist = r - b2.r;
        Type len = dist.Len();
        return dist * (-(Gamma * b2.m / (len * len * len)));
    }

    Body& operator=(Body body){
        m = body.m;
        r = body.r;
        v = body.v;
        a = body.a;
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

#endif //NBODIES_NBODIES_H
