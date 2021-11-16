#pragma once
#include <iostream>


class Boundary
{
public:
    Boundary(class Box* box_in);
    virtual void update() = 0;
    virtual bool correct_position() = 0;
    virtual bool correct_velocity() = 0;   // needed for reflective boundaries
    virtual bool correct_distance() = 0;   // needed for periodic boundaries
    //virtual double comp_volume() = 0;
    virtual ~Boundary() = default;

protected:
    class Box* box = nullptr;
};
