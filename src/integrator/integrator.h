#pragma once
#include <iostream>


class Integrator
{
public:
    Integrator(class Box *, double);
    virtual double next_step() = 0;
    virtual ~Integrator() = default;

    double dt;

protected:
    class Box* box = nullptr;
};
