#pragma once
#include <iostream>
#include <string>


class Boundary
{
public:
    Boundary(class Box *);
    virtual void update() = 0;
    virtual bool correct_position() = 0;
    virtual bool correct_velocity() = 0;   // needed for reflective boundaries
    virtual bool correct_distance() = 0;   // needed for periodic boundaries
    //virtual double comp_volume() = 0;
    virtual ~Boundary() = default;

    std::string label;

protected:
    class Box* box = nullptr;
};
