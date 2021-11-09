#pragma once
#include "integrator.h"


class VelocityVerlet : public Integrator
{
public:
    VelocityVerlet(class Box* box_in, double dt_in=0.01);
    void next_step();

private:
    double dt2;
    double ddt2;
};
