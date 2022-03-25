#pragma once
#include "integrator.h"


class VelocityVerlet : public Integrator
{
public:
    VelocityVerlet(class Box *, double = 0.01);
    double next_step();

private:
    double dt2;
    double ddt2;
};
