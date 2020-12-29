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

    void Update(Type step, vec<Type> da) {
        a += da;
        v += a * step;
        r += v * step;
    }

    double IteractSubtotalForce(const Body<Type> &b2) const{
        vec<Type> dist = b2.r - r;
        Type len = dist.Len();

        return Gamma * b2.m  / (len * len * len);
    }

    // For temporary utility bodies
    vec<Type> f(vec<Type> da, double step, double scale) {
        v += da * step * scale;
        r += v * step * scale;
        return a * r;
    }
};

#endif //NBODIES_NBODIES_H
