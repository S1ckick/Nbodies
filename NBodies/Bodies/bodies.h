//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_BODIES_H
#define NBODIES_BODIES_H

#include "../Math/vec.h"
#define G 6.6743e-11

template <class Type>
struct Body{
    Body(vec<Type> r, vec<Type> v,  Type m, vec<Type> a) : r(r), v(v), a(a), m(m) {};
    Body(vec<Type> r, vec<Type> v,  Type m) : r(r), v(v), m(m) {};
    Type m=0;
    vec<Type> r{0,0,0}, v{0,0,0}, a{0,0,0};
};

#endif //NBODIES_BODIES_H
