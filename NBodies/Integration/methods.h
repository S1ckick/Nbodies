//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_METHODS_H
#define NBODIES_METHODS_H

#include "../Nbodies/nbodies.h"

template <typename Type>
// b - {dist(b1, b2), v(v2), 0, |a| / dist}
vec<Type> RungeKutta4(const Body<Type> &b, double h) {
    vec<Type> k1, k2, k3, k4;
    Body<Type> tmp = b;
    k1 = b.a;
    k2 = tmp.f(k1, h, 0.5);
    k3 = tmp.f(k2, h, 0.5);
    k4 = tmp.f(k3, h, 1);

    return (k1 + k2 * 2 + k3 * 2 + k4) / 6;
}




#endif //NBODIES_METHODS_H
