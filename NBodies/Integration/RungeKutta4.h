//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_RUNGEKUTTA4_H
#define NBODIES_RUNGEKUTTA4_H


#include "../Bodies/bodies.h"

template<typename Type>
class RungeKutta4 {
public:
    vec<Type> SolveRungeKutta(Body<Type> external_body, Body<Type> target_body, Type f, double t);
    vec<Type> partial_step(vec<Type> f, vec<Type> df, Type scale, double t);
};

template <typename Type>
vec<Type> RungeKutta4<Type>::SolveRungeKutta(Body<Type> external_body, Body<Type> target_body, Type f, double t){
    vec<Type> k1{0,0,0}, k2{0,0,0}, k3{0,0,0}, k4{0,0,0}, velocity_update{0,0,0}, location_update{0,0,0};
//k1 - acceleration at current location
    k1 = (external_body.r - target_body.r) * f;

//k2 - acceleration 0.5 timesteps in the future based on k1 acceleration value
    velocity_update = partial_step(target_body.v, k1, 0.5, t);
    location_update = partial_step(target_body.r, velocity_update, 0.5, t);
    k2 = (external_body.r - location_update) * f;

//k3 acceleration 0.5 timesteps in the future using k2 acceleration
    velocity_update = partial_step(target_body.v, k2, 0.5, t);
    location_update = partial_step(target_body.r, velocity_update, 0.5, t);
    k3 = (external_body.r - location_update) * f;

//k4 - location 1 timestep in the future using k3 acceleration
    velocity_update = partial_step(target_body.v, k3, 1, t);
    location_update = partial_step(target_body.r, velocity_update, 1, t);
    k4 = (external_body.r - location_update) * f;

    return (k1 + k2 * 2 + k3 * 2 + k4) / 6;
}

template <typename Type>
vec<Type> RungeKutta4<Type>::partial_step(vec<Type> f, vec<Type> df, Type scale, double t){
    return f + df * t * scale;
}


#endif //NBODIES_RUNGEKUTTA4_H
