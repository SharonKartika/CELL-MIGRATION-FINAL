#pragma once
#include <random>

template <class T>
class Rand
{
public:
    static T getUnitRand()
    {
        std::random_device r;
        std::default_random_engine e1(r());
        std::uniform_real_distribution<T> unit_dist(0., 1.);
        T unitrand = unit_dist(e1);
        return unitrand;
    }
    static T getUniformRand(T a, T b)
    {
        std::random_device r;
        std::default_random_engine e1(r());
        std::uniform_real_distribution<T> unit_dist(a, b);
        return unit_dist(e1);
    }
};
