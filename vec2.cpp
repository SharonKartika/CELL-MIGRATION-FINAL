#include <cmath>

template <class T>
class VEC2
{
private:
public:
    T x;
    T y;
    VEC2()
    {
        this->x = 0;
        this->y = 0;
    }
    VEC2(T x, T y)
    {
        this->x = x;
        this->y = y;
    }
    VEC2 operator+(VEC2 const &obj)
    {
        VEC2<T> temp;
        temp.x = this->x + obj.x;
        temp.y = this->y + obj.y;
        return temp;
    }
    VEC2 operator-(VEC2 const &obj)
    {
        VEC2<T> temp;
        temp.x = this->x - obj.x;
        temp.y = this->y - obj.y;
        return temp;
    }
    VEC2 operator*(double const &a)
    {
        VEC2<T> temp;
        temp.x = x * a;
        temp.y = y * a;
        return temp;
    }
    void operator+=(VEC2 const &obj)
    {
        this->x += obj.x;
        this->y += obj.y;
    }
    void operator-=(VEC2 const &obj)
    {
        this->x -= obj.x;
        this->y -= obj.y;
    }
    VEC2 operator/(double const &a)
    {
        VEC2<T> temp;
        temp.x = this->x / a;
        temp.y = this->y / a;
        return temp;
    }
    double sqmag()
    {
        return pow(x, 2) + pow(y, 2);
    }
    double mag()
    {
        return sqrt(sqmag());
    }
    VEC2 unit()
    {
        return (*this) / mag();
    }
};
