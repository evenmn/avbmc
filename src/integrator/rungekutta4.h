#pragma once
#include "integrator.h"


class RungeKutta4 : public Integrator
{
public:
    RungeKutta4(class Box* box_in, double dt_in=0.01);
    void next_step();

private:
    mat pos_old, vel_old, pos_new, vel_new, acc_new;
    mat K1x, K1v, K2x, K2v, K3x, K3v, K4x, K4v;
    vec potengs;
};
