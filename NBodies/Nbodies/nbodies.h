
#ifndef NBODIES_NBODIES_H
#define NBODIES_NBODIES_H

#ifdef NUMBER_DOUBLE_DOUBLE
#include <qd/dd_real.h>
#endif

#include <random>

#include "../Utils/vec.h"
#include "summation.h"

template <typename Type>
static const Type Gamma(6.6743e-11);

template <class Type>
struct Body {
  Type m;
  vec<Type> r, v;

  Body(vec<Type> r1, vec<Type> v1, Type m1) : r(r1), v(v1), m(m1){};

  vec<Type> IteractSubtotalForce(const Body<Type> &b2) const {
    vec<Type> dist = r - b2.r;
    Type len = dist.Len();
    return dist * (-(b2.m / (len * len * len)));
  }

  Body &operator=(const Body &body) {
    r = body.r;
    v = body.v;
    m = body.m;
    return *this;
  }

  friend Body operator+(Body body1, const Body &body2) {
    body1.r = body1.r + body2.r;
    body1.v = body1.v + body2.v;
    return body1;
  }

  friend Body operator*(Body body1, Type number) {
    body1.r = body1.r * number;
    body1.v = body1.v * number;
    return body1;
  }
};

template <typename Type>
vec<Type> force(const vec<Type> &v1, const vec<Type> &v2, Type m1, Type m2) {
  vec<Type> dist = v1 - v2;
  Type len = dist.Len();
  return dist * (-(m1 * m2 / (Gamma<Type> * len * len * len)));
}

template <typename Type>
void copyBodies(const std::vector<Body<Type>> &bodies,
                std::vector<Body<Type>> &newBodies) {
  newBodies.clear();
  newBodies.resize(bodies.size());
  for (int i = 0; i < bodies.size(); i++) {
    newBodies.push_back(bodies[i]);
  }
}

template <typename Type>
std::vector<Body<Type>> operator+(const std::vector<Body<Type>> &bodies_1,
                                  const std::vector<Body<Type>> &bodies_2) {
  std::vector<Body<Type>> res;
  res.reserve(bodies_1.size());
  for (int i = 0; i < bodies_1.size(); i++) {
    Body<Type> item = bodies_1[i] + bodies_2[i];
    res.push_back(item);
  }
  return res;
}

template <typename Type>
std::vector<Body<Type>> operator*(const std::vector<Body<Type>> &bodies_1,
                                  Type number) {
  std::vector<Body<Type>> res;
  res.reserve(bodies_1.size());
  for (int i = 0; i < bodies_1.size(); i++) {
    Body<Type> item = bodies_1[i] * number;
    res.push_back(item);
  }
  return res;
}

template <typename Type>
struct acceleration_proxy {
  const std::vector<Body<Type>> &bodies;
  const size_t bodyind;
  explicit acceleration_proxy(const std::vector<Body<Type>> &bodies, size_t bi)
      : bodies(bodies), bodyind(bi) {}
  vec<Type> operator[](size_t n) const {
    if (n == bodyind) return vec<Type>(0);
    return bodies[bodyind].IteractSubtotalForce(bodies[n]);
  }
};

template <typename Type>
struct total_mass_proxy {
  const std::vector<Body<Type>> &bodies;
  explicit total_mass_proxy(const std::vector<Body<Type>> &bodies)
      : bodies(bodies) {}
  Type operator[](size_t n) const { return bodies[n].m; }
};

// derivative
template <typename Type>
void f(const std::vector<Body<Type>> &bodies,
       std::vector<Body<Type>> &fbodies) {
  for (int body_1_idx = 0; body_1_idx < bodies.size(); body_1_idx++) {
    acceleration_proxy<Type> a(bodies, body_1_idx);
    vec<Type> total_force =
        summation<vec<Type>, acceleration_proxy<Type>>(a, bodies.size());
    fbodies[body_1_idx].r = bodies[body_1_idx].v;
    fbodies[body_1_idx].v = total_force;
  }
}

#endif  // NBODIES_NBODIES_H
