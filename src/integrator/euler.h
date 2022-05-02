#pragma once
#include "integrator.h"


class Euler : public Integrator
{
public:
    Euler(class Box *, double = 0.001);
    double next_step() override;
};
