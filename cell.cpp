#include "vec2.cpp"
#include "random.cpp"
#include "math.h"
#include <iostream>
template <class T>
class Cell
{
private:
public:
    VEC2<T> pos, vel, acc, eta;
    Cell()
    {
        pos = VEC2<T>(10, -10);
        vel = VEC2<T>(0, 0);
        // vel = VEC2<T>(Rand<T>::getUniformRand(-4, 4), Rand<T>::getUniformRand(-10, 10));
        acc = VEC2<T>(0, 0);
        T theta = Rand<T>::getUniformRand(0, 2*M_PI);
        eta = VEC2(cos(theta), sin(theta));
    }
    void showPos()
    {
        std::cout << pos.x << " " << pos.y << std::endl;
    }
    void setRandomPosition(T width, T height)
    {
        T x = Rand<T>::getUniformRand(-width, width);
        T y = Rand<T>::getUniformRand(-height, height);
        pos.x = x;
        pos.y = y;
    }
    void update(T dt)
    {
        vel += acc * dt;
        pos += vel * dt;
        acc = VEC2<T>(0, 0);
    }
};