//
// Created by Максим on 01.01.2021.
//

#ifndef NBODIES_SUMMATION_H
#define NBODIES_SUMMATION_H

#include "energy.h"
//#include <qd/dd_real.h>



template<class Type, class A>
Type summation_neumaier(const A &container, size_t begin, size_t end, Type& correction){
    Type s = Type(0), c = Type(0), t = Type(0);
    for(size_t i = begin; i < end; ++i){
        Type item = container[i];
        t = s + item;
        if(abs(s) >= abs(item)) {
            c += (s - t) + item;
        }
        else {
            c += (item - t) + s;
        }
        s = t;
    }
    correction = c;
    return s + correction;
}

template <class Type>
Type summation_kahan(const Type &a, const Type &b, Type &correction){
    Type corrected = b - correction;
    Type new_sum = a + corrected;
    correction = (new_sum - a) - corrected;
    return new_sum;
}

template <class Type, class A>
Type summation_arr(A container, size_t begin, size_t end, Type* correction){
    if(begin == end){
        return Type(0);
    }
    Type sum;
    sum = container[begin];
    for(size_t i = begin + 1; i < end; ++i){
        sum = summation_kahan(sum, container[i], *correction);
    }
    return sum;
}

template <class Type, class A>
Type summation(const A& container, size_t size){
    Type correction(0);

    return summation_neumaier(container, 0, size, correction);
}

#endif //NBODIES_SUMMATION_H
