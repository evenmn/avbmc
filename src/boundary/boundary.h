#pragma once
#include <iostream>
#include <string>
#include <valarray>


class Boundary
{
public:
    Boundary(class Box *);
    virtual inline void correct_position(unsigned int) {};
    virtual inline void correct_velocity(unsigned int) {};              // needed for reflective boundaries
    virtual inline void correct_distance(std::valarray<double> &) {};   // needed for periodic boundaries
    //virtual double comp_volume() = 0;
    virtual ~Boundary() = default;

    std::string label;

protected:
    class Box* box = nullptr;
};
