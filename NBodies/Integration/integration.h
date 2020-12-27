//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_INTEGRATION_H
#define NBODIES_INTEGRATION_H

#include "../Bodies/bodies.h"
#include "RungeKutta4.h"

#define G 6.6743e-11
#include <vector>

template <typename Type>
class Integration{
public:
    void compute(std::vector<Body<Type>> &bodies, double m_time_step);
};

template <typename Type>
void Integration<Type>::compute(std::vector< Body<Type> > &bodies, double m_time_step){
    RungeKutta4<Type> rk;
    for(int i = 0; i < bodies.size(); i++){
        Type mass1 = bodies[i].m;
        vec<Type> acceleration;
        for(int j = i + 1; j < bodies.size(); j++){
            Type mass2 = bodies[j].m;

            vec<Type> dist = bodies[i].r - bodies[j].r;

            Type f = G * (mass1 * mass2) / (dist.Len2());

            acceleration+=rk.SolveRungeKutta(bodies[j], bodies[i], f, m_time_step);

        }
        bodies[i].v += acceleration * m_time_step;

    }
    for(int i = 0; i < bodies.size(); i++){
        bodies[i].r += bodies[i].v * m_time_step;
    }
}




#endif //NBODIES_INTEGRATION_H
