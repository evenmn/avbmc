#pragma once
#include <iostream>
#include <string>


class Boundary
{
public:
    Boundary(class Box *);
    virtual bool correct_position();
    bool correct_velocity();   // needed for reflective boundaries
    bool correct_distance();   // needed for periodic boundaries
    //virtual double comp_volume() = 0;
    virtual ~Boundary() = default;

    std::string label;

protected:
    class Box* box = nullptr;
};
