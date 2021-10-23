#pragma once
#include "integrator.h"


class Euler : public Integrator
{
public:
    Euler(class Box* box_in, double dt_in=0.001);
    void next_step();
};
