//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_VEC_H
#define NBODIES_VEC_H

#ifdef NUMBER_DOUBLE_DOUBLE
#include <qd/dd_real.h>
#endif

#ifdef NUMBER_DOUBLE
#include <cmath>
#endif


template <class type>
class vec {
public:
    type X, Y, Z;

    vec() : X(type(0)), Y(type(0)), Z(type(0)) {}
    vec(type A, type B, type C) : X(A), Y(B), Z(C) {}
    vec(const vec& v) : X(v.X), Y(v.Y), Z(v.Z) {}
    explicit vec(type A) : X(A), Y(A), Z(A) {}

    operator type *(void) { return &X; }

    inline bool operator==(const vec &V) const {
        return (X == V.X && Y == V.Y && Z == V.Z);
    }

    inline bool operator!=(const vec &V) const {
        return !(X == V.X && Y == V.Y && Z == V.Z);
    }

    inline vec operator+(const vec &V) const {
        return vec(X + V.X, Y + V.Y, Z + V.Z);
    }

    inline vec operator+(type N) const {
        return vec(X + N, Y + N, Z + N);
    }

    inline vec &operator+=(const vec &V) {
        X += V.X;
        Y += V.Y;
        Z += V.Z;
        return *this;
    }

    inline vec &operator+=(type N) {
        X += N;
        Y += N;
        Z += N;
        return *this;
    }

    inline vec operator-(const vec &V) const {
        return vec(X - V.X, Y - V.Y, Z - V.Z);
    }

    inline vec operator-(type N) const {
        return vec(X - N, Y - N, Z - N);
    }

    inline vec &operator-=(const vec &V) {
        X -= V.X;
        Y -= V.Y;
        Z -= V.Z;
        return *this;
    }

    inline vec &operator-=(type N) {
        X -= N;
        Y -= N;
        Z -= N;
        return *this;
    }

    inline vec operator/(const vec &V) const {
        if (V.X == 0 || V.Y == 0 || V.Z == 0) return vec(X, Y, Z);
        return vec(X / V.X, Y / V.Y, Z / V.Z);
    }

    inline vec operator/(type N) const {
        if (N == 0) return vec(X, Y, X);
        return vec(X / N, Y / N, Z / N);
    }

    inline vec &operator/=(const vec &V) {
        if (V.X == 0 || V.Y == 0 || V.Z == 0) return *this;
        X /= V.X;
        Y /= V.Y;
        Z /= V.Z;
        return *this;
    }

    inline vec &operator/=(type N) {
        if (N == 0) return *this;
        X /= N;
        Y /= N;
        Z /= N;
        return *this;
    }

    inline vec operator*(const vec &N) const {
        return vec(X * N.X, Y * N.Y, Z * N.Z);
    }

    inline vec operator*(type N) const {
        return vec(X * N, Y * N, Z * N);
    }

    inline vec &operator*=(const vec &V) {
        X *= V.X;
        Y *= V.Y;
        Z *= V.Z;
        return *this;
    }

    inline vec &operator*=(type N) {
        X *= N;
        Y *= N;
        Z *= N;
        return *this;
    }

    inline vec operator-() const {
        return vec(-X, -Y, -Z);
    }

    inline vec &operator-() {
        X = -X;
        Y = -Y;
        Z = -Z;
        return *this;
    }

    inline type operator&(const vec &V) const {
        return X * V.X + Y * V.Y + Z * V.Z;
    }

    inline vec operator%(const vec &V) const {
        return vec(Y * V.Z - Z * V.Y, Z * V.X - X * V.Z, X * V.Y - Y * V.X);
    }

    inline vec operator%=(const vec &V) {
        X = Y * V.Z - V.Y * Z;
        Y = V.X * Z - X * V.Z;
        Z = X * V.Y - V.X * Y;
        return *this;
    }

    inline vec Normalizing() const {
        type len = *this & *this;
        if (len != 0 && len != 1) {
            len = sqrt(len);
            return vec(X / len, Y / len, Z / len);
        }
        return *this;
    }

    inline vec &Normalize() {
        type len = *this & *this;
        if (len != 0 && len != 1) {
            len = sqrt(len);
            X /= len;
            Y /= len;
            Z /= len;
        }
        return *this;
    }

    inline type Len() const {
        return sqrt(X*X + Y*Y + Z*Z);
    }

    inline type Len2() const {
        return type(X * X + Y * Y + Z * Z);
    }

    inline type &operator[](int I) {
        if (I >= 0 && I < 3)
            return *(&X + I);
        return 0;
    }
};

template <class Type>
Type abs(const vec<Type> &v) {
  return v.Len2();
}

#endif //NBODIES_VEC_H
