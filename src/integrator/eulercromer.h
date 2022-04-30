#pragma once
#include "integrator.h"


class EulerCromer : public Integrator
{
public:
    EulerCromer(class Box *, double = 0.01);
    double next_step() override;
};
