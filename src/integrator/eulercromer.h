#pragma once
#include "integrator.h"


class EulerCromer : public Integrator
{
public:
    EulerCromer(class Box* box_in, double dt_in=0.01);
    void next_step();
};
